#include <stdint.h>

static int64_t rule_calls;
static int counting_enabled = 1;

void fortnum_enzyme_rule_counter_record(void)
{
    if (counting_enabled) {
        ++rule_calls;
    }
}

void fortnum_enzyme_rule_counter_reset(void)
{
    rule_calls = 0;
    counting_enabled = 1;
}

int64_t fortnum_enzyme_rule_counter_calls(void)
{
    return rule_calls;
}

void fortnum_enzyme_rule_counter_disable(void)
{
    counting_enabled = 0;
}
