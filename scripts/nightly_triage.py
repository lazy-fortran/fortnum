#!/usr/bin/env python3
"""Run the nightly jobs and triage their failures.

Subcommands, all run inside GitHub Actions:

  plan            turn the dispatch input into the list of jobs to run
  run JOB         run one nightly job's command
  bisect JOB      after JOB failed: bisect its command between the last
                  commit where it passed and this commit, write the result
  report          file or update the nightly-failure issue, or close it
                  once every job it records is green again

`bisect` runs inside the failed job, so the job's toolchain is already set
up; only the command below is repeated at each commit. A nonzero exit,
including a failed build, counts as bad.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# key: job name as shown in nightly.yml, the job name it had in the merge
# queue (to find a passing commit before the first green nightly), and the
# command repeated at each bisection step.
JOBS = {
    "enzyme": {
        "name": "released Flang 22 + pinned Enzyme",
        "legacy": ("ci.yml", "released Flang 22 + pinned Enzyme"),
        "command": r"""
wrappers="$RUNNER_TEMP/enzyme-wrappers"
rm -rf "$wrappers" && mkdir -p "$wrappers"
if [ ! -d ../fortsym/.git ]; then
    git clone https://github.com/lazy-fortran/fortsym.git ../fortsym
fi
# fpm builds every generator in tools/codegen together, so the tree
# must satisfy every lock, even though only the Enzyme wrappers are used.
tools/codegen/checkout_fortsym.sh ../fortsym
(
    cd tools/codegen
    FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR="$wrappers" \
        fo exec gen_enzyme_scalar_wrappers
)
cmake -S . -B build-enzyme -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFORTNUM_BUILD_TESTING=ON \
    -DFORTNUM_BUILD_EXAMPLES=OFF \
    -DFORTNUM_ENABLE_ENZYME=ON \
    -DFORTNUM_ENZYME_REQUIRED=ON \
    -DFORTNUM_ENZYME_PLUGIN="$FORTNUM_CI_ENZYME_PLUGIN" \
    -DFORTNUM_ENZYME_WRAPPER_DIR="$wrappers"
cmake --build build-enzyme -j 2
ctest --test-dir build-enzyme --output-on-failure -R enzyme
""",
    },
}

LABEL = "nightly-failure"
MARKER = re.compile(r"<!-- nightly-failing-jobs: ([a-z0-9_,-]*) -->")
FAILED = ("failure", "timed_out")


def parse_jobs(text):
    """Job keys from a dispatch input; empty means every job."""
    text = (text or "").strip()
    if not text or text == "all":
        return list(JOBS)
    keys = []
    for item in text.split(","):
        item = item.strip()
        if item and item not in keys:
            if item not in JOBS:
                raise ValueError(f"unknown job {item!r}; known: {', '.join(JOBS)}")
            keys.append(item)
    if not keys:
        raise ValueError(f"no job in {text!r}")
    return keys


def failing_jobs(jobs):
    """Keys of the nightly jobs that failed, from the run's job list."""
    by_name = {spec["name"]: key for key, spec in JOBS.items()}
    return sorted(by_name[j["name"]] for j in jobs
                  if j["name"] in by_name and j.get("conclusion") in FAILED)


def job_passed(jobs, name):
    return any(j["name"] == name and j.get("conclusion") == "success" for j in jobs)


def issue_title(branch):
    return f"Nightly failure on {branch}"


def recorded_jobs(texts):
    """All failing jobs recorded in an issue body and its comments."""
    keys = set()
    for text in texts:
        for match in MARKER.finditer(text or ""):
            keys.update(k for k in match.group(1).split(",") if k)
    return keys


def should_close(green, recorded):
    return bool(recorded) and set(recorded) <= set(green)


def plain(text):
    """Commit subjects go in a code span, so they cannot mention anyone."""
    return "`" + text.replace("`", "'") + "`"


def failure_report(*, branch, bad, run_url, server, repo, failures):
    """Issue text; failures maps a job key to its URL and bisection result."""
    lines = [f"The nightly run failed on `{branch}` at {bad}.",
             f"Run: {run_url}", ""]
    for key in sorted(failures):
        item = failures[key]
        result = item.get("result") or {}
        lines.append(f"- **{JOBS[key]['name']}** ([log]({item['url']}))")
        good = result.get("good")
        if good:
            lines.append(f"  - Commit range: {good}..{bad} "
                         f"({result.get('commits', '?')} commits, "
                         f"[compare]({server}/{repo}/compare/{good}...{bad}))")
        elif result:
            lines.append("  - Commit range: unknown, no earlier passing run of this job")
        else:
            lines.append("  - Not bisected: the job failed before its command ran, "
                         "or the bisection did not finish")
        culprit = result.get("culprit")
        if culprit:
            via = f" from #{item['pull']}" if item.get("pull") else ""
            lines.append(f"  - First bad commit: {culprit['sha']} "
                         f"{plain(culprit['subject'])}{via}")
        if result.get("note"):
            lines.append(f"  - {result['note']}")
    lines.append("")
    lines.append(f"<!-- nightly-failing-jobs: {','.join(sorted(failures))} -->")
    return "\n".join(lines)


# ---- GitHub and git plumbing -------------------------------------------------

def env(name):
    return os.environ[name]


def api(path, method="GET", data=None):
    url = path if path.startswith("http") else f"{env('GITHUB_API_URL')}/{path}"
    body = json.dumps(data).encode() if data is not None else None
    request = urllib.request.Request(url, data=body, method=method, headers={
        "Authorization": f"Bearer {env('GH_TOKEN')}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28"})
    with urllib.request.urlopen(request) as response:
        text = response.read().decode()
    return json.loads(text) if text.strip() else None


def run_jobs(repo, run_id):
    return api(f"repos/{repo}/actions/runs/{run_id}/jobs?per_page=100")["jobs"]


def git(*args, check=True):
    return subprocess.run(["git", *args], check=check, capture_output=True,
                          text=True).stdout.strip()


def is_ancestor(old, new):
    return subprocess.run(["git", "merge-base", "--is-ancestor", old, new],
                          capture_output=True).returncode == 0


def last_green(repo, branch, key, current_run, bad):
    """Newest commit, an ancestor of bad, at which this job passed.

    Nightly runs on the branch come first. Before the first green nightly,
    the job's earlier runs in the merge-queue workflow count.
    """
    spec = JOBS[key]
    query = urllib.parse.urlencode({"branch": branch, "per_page": 30})
    runs = api(f"repos/{repo}/actions/workflows/nightly.yml/runs?{query}")
    for run in runs["workflow_runs"]:
        if run["id"] == current_run or run["status"] != "completed":
            continue
        if (job_passed(run_jobs(repo, run["id"]), spec["name"])
                and is_ancestor(run["head_sha"], bad)):
            return run["head_sha"]
    workflow, name = spec["legacy"]
    runs = api(f"repos/{repo}/actions/workflows/{workflow}/runs?per_page=30")
    for run in runs["workflow_runs"]:
        if run["event"] not in ("merge_group", "push", "workflow_dispatch"):
            continue
        if run["status"] != "completed" or not is_ancestor(run["head_sha"], bad):
            continue
        if job_passed(run_jobs(repo, run["id"]), name):
            return run["head_sha"]
    return None


def run_command(key):
    status = subprocess.run(["bash", "-euo", "pipefail", "-c",
                             JOBS[key]["command"]]).returncode
    return 0 if status == 0 else 1


def cmd_plan(args):
    keys = parse_jobs(args.jobs)
    with open(env("GITHUB_OUTPUT"), "a") as out:
        out.write(f"jobs={json.dumps(keys)}\n")


def cmd_run(args):
    return run_command(args.job)


def cmd_step(args):
    """One `git bisect run` step: 0 good, 1 bad."""
    if Path(".gitmodules").exists():
        git("submodule", "update", "--init", "--recursive", check=False)
    status = run_command(args.job)
    git("reset", "--hard", "--quiet", check=False)  # generators may edit tracked files
    return status


def bisect(key, good, bad):
    """First bad commit in good..bad, or None if the job fails at good too."""
    script = Path(os.environ.get("RUNNER_TEMP", "/tmp")) / "nightly_triage.py"
    shutil.copy(__file__, script)  # bisection checks out older versions of this file
    step = [sys.executable, str(script), "step", key]
    git("checkout", "--quiet", "--force", good)
    try:
        good_still_passes = subprocess.run(step).returncode == 0
    finally:
        git("checkout", "--quiet", "--force", bad)
    if not good_still_passes:
        return None
    git("bisect", "start", bad, good)
    try:
        subprocess.run(["git", "bisect", "run", *step], check=True)
        sha = git("rev-parse", "refs/bisect/bad")
    finally:
        print(git("bisect", "log", check=False))
        git("bisect", "reset", check=False)
    return {"sha": sha, "subject": git("log", "-1", "--format=%s", sha)}


def cmd_bisect(args):
    repo, run_id = env("GITHUB_REPOSITORY"), int(env("GITHUB_RUN_ID"))
    branch, bad = env("GITHUB_REF_NAME"), env("GITHUB_SHA")
    good = last_green(repo, branch, args.job, run_id, bad)
    result = {"good": good, "commits": 0, "culprit": None, "note": None}
    if not good:
        pass
    elif good == bad:
        result["note"] = ("No commit since the job last passed: "
                          "the failure is not from a code change.")
    else:
        result["commits"] = int(git("rev-list", "--count", f"{good}..{bad}"))
        result["culprit"] = bisect(args.job, good, bad)
        if not result["culprit"]:
            result["note"] = ("The job also fails at the last passing commit when "
                              "rerun, so the failure is not from a commit in the range.")
    out = Path(args.results)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{args.job}.json").write_text(json.dumps(result))
    print(json.dumps(result, indent=2))


def open_issue(repo, branch):
    query = urllib.parse.urlencode({"labels": LABEL, "state": "open", "per_page": 50})
    for issue in api(f"repos/{repo}/issues?{query}"):
        if issue["title"] == issue_title(branch) and "pull_request" not in issue:
            return issue
    return None


def ensure_label(repo):
    try:
        api(f"repos/{repo}/labels/{LABEL}")
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        api(f"repos/{repo}/labels", "POST", {
            "name": LABEL, "color": "B60205",
            "description": "A nightly job failed on the default branch"})


def cmd_report(args):
    repo, run_id = env("GITHUB_REPOSITORY"), int(env("GITHUB_RUN_ID"))
    server, branch, bad = env("GITHUB_SERVER_URL"), env("GITHUB_REF_NAME"), env("GITHUB_SHA")
    run_url = f"{server}/{repo}/actions/runs/{run_id}"
    jobs = run_jobs(repo, run_id)
    failing = failing_jobs(jobs)
    issue = open_issue(repo, branch)
    if not failing:
        planned = parse_jobs(args.jobs)
        green = [k for k in planned if job_passed(jobs, JOBS[k]["name"])]
        if not issue:
            return 0
        comments = api(issue["comments_url"] + "?per_page=100")
        recorded = recorded_jobs([issue["body"]] + [c["body"] for c in comments])
        if not should_close(green, recorded):
            print(f"issue #{issue['number']} records {sorted(recorded)}; green here: {green}")
            return 0
        api(f"repos/{repo}/issues/{issue['number']}/comments", "POST",
            {"body": f"Every recorded nightly job passes again at {bad}: {run_url}"})
        api(f"repos/{repo}/issues/{issue['number']}", "PATCH",
            {"state": "closed", "state_reason": "completed"})
        return 0
    failures = {}
    for key in failing:
        item = {"url": next(j["html_url"] for j in jobs if j["name"] == JOBS[key]["name"])}
        path = Path(args.results) / f"{key}.json"
        if path.exists():
            item["result"] = json.loads(path.read_text())
            culprit = item["result"].get("culprit")
            if culprit:
                pulls = api(f"repos/{repo}/commits/{culprit['sha']}/pulls")
                item["pull"] = pulls[0]["number"] if pulls else None
        failures[key] = item
    body = failure_report(branch=branch, bad=bad, run_url=run_url, server=server,
                          repo=repo, failures=failures)
    print(body)
    if issue:
        api(f"repos/{repo}/issues/{issue['number']}/comments", "POST", {"body": body})
    else:
        ensure_label(repo)
        api(f"repos/{repo}/issues", "POST",
            {"title": issue_title(branch), "body": body, "labels": [LABEL]})
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["plan", "run", "step", "bisect", "report"])
    parser.add_argument("job", nargs="?", choices=list(JOBS))
    parser.add_argument("--jobs", default="", help="dispatch input for plan and report")
    parser.add_argument("--results", default=os.path.join(
        os.environ.get("RUNNER_TEMP", "/tmp"), "nightly-triage"))
    args = parser.parse_args(argv)
    if args.action in ("run", "step", "bisect") and not args.job:
        parser.error(f"{args.action} needs a job")
    return {"plan": cmd_plan, "run": cmd_run, "step": cmd_step, "bisect": cmd_bisect,
            "report": cmd_report}[args.action](args) or 0


if __name__ == "__main__":
    sys.exit(main())
