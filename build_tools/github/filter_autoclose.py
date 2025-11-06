import os
from datetime import datetime, timezone

import requests

from build_tools.get_comment import get_headers

CUTOFF_DAYS = 14


def _get_paginated_results(url, headers):
    results = []
    page = 1
    while True:
        paged_url = f"{url}?per_page=100&page={page}"
        response = requests.get(paged_url, headers=headers)
        response.raise_for_status()
        data = response.json()
        if not data:
            break
        results.extend(data)
        page += 1
    return results


def get_pull_requests_to_autoclose(GH_REPO, GITHUB_TOKEN, ALL_LABELED_PRS):
    # get "autoclose" labeled PRs that are older than 14 days
    all_labeled_prs = [int(x) for x in ALL_LABELED_PRS.split()]
    now = datetime.now(timezone.utc)
    pull_requests_to_autoclose = []

    for pr in all_labeled_prs:
        timestamp_label_last_set = None
        url = f"https://api.github.com/repos/{GH_REPO}/issues/{pr}/events"
        events = _get_paginated_results(url, get_headers(GITHUB_TOKEN))
        for event in events:
            if event["event"] == "labeled" and event["label"]["name"] == "autoclose":
                timestamp_label_last_set = event["created_at"]
                timestamp_label_last_set = datetime.strptime(
                    timestamp_label_last_set, "%Y-%m-%dT%H:%M:%SZ"
                )
                timestamp_label_last_set = timestamp_label_last_set.replace(
                    tzinfo=timezone.utc
                )
        if timestamp_label_last_set is None:
            # we only filter through PRs that have an autoclose label set, so a failure
            # to find that key would mean something is off and the action should fail:
            raise KeyError("Could not find 'autoclose' label in pre-filtered PRs.")
        if (now - timestamp_label_last_set).days > 14:
            pull_requests_to_autoclose.append(pr)

    return ",".join([str(x) for x in pull_requests_to_autoclose])


if __name__ == "__main__":
    GH_REPO = os.getenv("GH_REPO", "")
    GITHUB_TOKEN = os.getenv("GITHUB_TOKEN", "")
    ALL_LABELED_PRS = os.getenv("ALL_LABELED_PRS", "")

    pull_requests_to_autoclose = get_pull_requests_to_autoclose(
        GH_REPO, GITHUB_TOKEN, ALL_LABELED_PRS
    )

    with open(os.getenv("GITHUB_ENV"), "a") as f:
        f.write(f"PULL_REQUESTS={pull_requests_to_autoclose}\n")
