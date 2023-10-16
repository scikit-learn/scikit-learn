# Update tracking issue if Cirrus fails nightly job

# if [[ "$CIRRUS_CRON" != "nightly" ]]; then
#     exit 0
# fi

# TEST_PASSED is either "true" or "false"
TEST_PASSED="$1"

python -m venv .venv
source .venv/bin/activate
python -m pip install defusedxml PyGithub

LINK_TO_RUN="https://cirrus-ci.com/build/$CIRRUS_BUILD_ID"

echo $CIRRUS_TASK_NAME
echo $CIRRUS_REPO_FULL_NAME
echo $LINK_TO_RUN
echo $TEST_PASSED

python maint_tools/update_tracking_issue.py \
    $BOT_GITHUB_TOKEN \
    $CIRRUS_TASK_NAME \
    $CIRRUS_REPO_FULL_NAME \
    $LINK_TO_RUN \
    --tests-passed $TEST_PASSED \
    --auto-close false
