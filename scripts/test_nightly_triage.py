"""Guards for the nightly triage logic, with no GitHub access."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nightly_triage as nt  # noqa: E402

KEY = next(iter(nt.JOBS))
NAME = nt.JOBS[KEY]["name"]


class TriageTests(unittest.TestCase):
    def test_job_selection(self):
        self.assertEqual(nt.parse_jobs(""), list(nt.JOBS))
        self.assertEqual(nt.parse_jobs(f" {KEY},{KEY}"), [KEY])
        with self.assertRaises(ValueError):
            nt.parse_jobs("no-such-job")

    def test_failing_jobs_from_run(self):
        jobs = [{"name": NAME, "conclusion": "failure"},
                {"name": "Plan nightly jobs", "conclusion": "failure"}]
        self.assertEqual(nt.failing_jobs(jobs), [KEY])
        jobs[0]["conclusion"] = "cancelled"
        self.assertEqual(nt.failing_jobs(jobs), [])
        self.assertFalse(nt.job_passed(jobs, NAME))
        jobs[0]["conclusion"] = "success"
        self.assertTrue(nt.job_passed(jobs, NAME))

    def test_report_and_closing(self):
        result = {"good": "a" * 40, "commits": 3, "note": None,
                  "culprit": {"sha": "c" * 40, "subject": "Break it, thanks @someone"}}
        body = nt.failure_report(
            branch="main", bad="b" * 40, run_url="https://x/run",
            server="https://github.com", repo="o/r",
            failures={KEY: {"url": "https://x/job", "result": result, "pull": 42}})
        self.assertIn(f"**{NAME}** ([log](https://x/job))", body)
        self.assertIn("compare/" + "a" * 40 + "..." + "b" * 40, body)
        self.assertIn("`Break it, thanks @someone` from #42", body)
        outside_code = body.replace("`Break it, thanks @someone`", "")
        self.assertNotIn("@", outside_code)
        recorded = nt.recorded_jobs([body, "a comment"])
        self.assertEqual(recorded, {KEY})
        self.assertFalse(nt.should_close([], recorded))
        self.assertTrue(nt.should_close(list(nt.JOBS), recorded))
        self.assertFalse(nt.should_close(list(nt.JOBS), set()))

    def test_bisect_finds_the_breaking_commit(self):
        """git bisect over a real repository, with a probe job added."""
        with tempfile.TemporaryDirectory() as tmp:
            repo, tools = Path(tmp, "repo"), Path(tmp, "tools")
            repo.mkdir()
            tools.mkdir()

            def git(*args):
                return subprocess.run(["git", *args], cwd=repo, check=True,
                                      capture_output=True, text=True).stdout.strip()
            git("init", "-q", "-b", "main")
            git("config", "user.email", "probe@example.org")
            git("config", "user.name", "probe")
            git("config", "commit.gpgsign", "false")
            git("config", "core.hooksPath", str(tools))  # an empty hook directory
            shas = []
            for n, state in enumerate(["ok", "ok", "ok", "broken", "broken"]):
                Path(repo, "state").write_text(state)
                git("add", "state")
                git("commit", "-q", "--allow-empty", "-m", f"commit {n}")
                shas.append(git("rev-parse", "HEAD"))
            text = Path(nt.__file__).read_text().replace(
                "JOBS = {", 'JOBS = {"probe": {"name": "probe", "legacy": ("ci.yml", '
                '"probe"), "command": "grep -qx ok state"},', 1)
            Path(tools, "nightly_triage.py").write_text(text)

            def bisect(good):
                code = ("import json, sys; sys.path.insert(0, sys.argv[1]); "
                        "import nightly_triage as t; "
                        "print(json.dumps(t.bisect('probe', sys.argv[2], sys.argv[3])))")
                out = subprocess.run(
                    [sys.executable, "-c", code, str(tools), good, shas[-1]],
                    cwd=repo, check=True, capture_output=True, text=True,
                    env={"PATH": "/usr/bin:/bin", "RUNNER_TEMP": tmp})
                return json.loads(out.stdout.strip().splitlines()[-1])

            culprit = bisect(shas[0])
            self.assertEqual(culprit["sha"], shas[3])
            self.assertEqual(culprit["subject"], "commit 3")
            self.assertEqual(git("rev-parse", "HEAD"), shas[-1])
            # A job that fails at the passing commit as well has no culprit.
            self.assertIsNone(bisect(shas[3]))


if __name__ == "__main__":
    unittest.main()
