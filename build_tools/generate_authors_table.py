"""
This script generates an HTML table of contributors, with names and avatars.
The list is generated from scikit-learn's teams on GitHub, plus a small number
of hard-coded contributors.

The table should be updated for each new inclusion in the teams.
Generating the table requires admin rights.
"""

import getpass
import sys
import time
from pathlib import Path
import requests

# Prompt for GitHub credentials
print("Input user:", file=sys.stderr)
user = input()
token = getpass.getpass("Input access token:\n")
auth = (user, token)

# Create a session to reuse connections and set authentication globally.
session = requests.Session()
session.auth = auth

LOGO_URL = "https://avatars2.githubusercontent.com/u/365630?v=4"
REPO_FOLDER = Path(__file__).resolve().parent.parent


def get(url):
    """Wrapper for GET requests with rate limit handling."""
    # Use increasing sleep times for rate-limit retries.
    for sleep_time in [10, 30, 60]:
        reply = session.get(url)
        try:
            json_resp = reply.json()
        except ValueError:
            json_resp = {}
        # Check if the API rate limit is exceeded.
        if not ("message" in json_resp and "API rate limit exceeded" in json_resp.get("message", "")):
            break
        print("API rate limit exceeded, waiting...")
        time.sleep(sleep_time)

    reply.raise_for_status()
    return reply


def get_contributors():
    """Get the list of contributor profiles. Requires admin rights."""
    # Define team lists.
    core_devs = []
    documentation_team = []
    contributor_experience_team = []
    comm_team = []
    # Define team slugs.
    core_devs_slug = "core-devs"
    contributor_experience_team_slug = "contributor-experience-team"
    comm_team_slug = "communication-team"
    documentation_team_slug = "documentation-team"

    entry_point = "https://api.github.com/orgs/scikit-learn/"

    # Retrieve team members.
    for team_slug, lst in zip(
        (core_devs_slug, contributor_experience_team_slug, comm_team_slug, documentation_team_slug),
        (core_devs, contributor_experience_team, comm_team, documentation_team),
    ):
        print(f"Retrieving {team_slug}")
        for page in range(1, 3):  # Pages 1 and 2 (30 per page)
            reply = get(f"{entry_point}teams/{team_slug}/members?page={page}")
            lst.extend(reply.json())

    # Retrieve all members.
    print("Retrieving members")
    members = []
    for page in range(1, 4):  # Pages 1, 2, and 3
        reply = get(f"{entry_point}members?page={page}")
        members.extend(reply.json())

    # Keep only the login names.
    core_devs = {c["login"] for c in core_devs}
    documentation_team = {c["login"] for c in documentation_team}
    contributor_experience_team = {c["login"] for c in contributor_experience_team}
    comm_team = {c["login"] for c in comm_team}
    members = {c["login"] for c in members}

    # Add missing contributors.
    members |= {"dubourg", "mbrucher", "thouis", "jarrodmillman"}
    members |= {"Angel Soler Gollonet"}
    # Remove CI bots.
    members -= {"sklearn-ci", "sklearn-wheels", "sklearn-lgtm"}
    contributor_experience_team -= core_devs  # Remove ogrisel from contributor_experience_team

    # Determine emeritus contributors.
    emeritus = members - core_devs - contributor_experience_team - comm_team - documentation_team

    # Hard-coded emeritus teams.
    emeritus_contributor_experience_team = {"cmarmo"}
    emeritus_comm_team = {"reshamas"}

    # Exclude hard-coded emeritus teams from the original emeritus.
    emeritus -= emeritus_contributor_experience_team | emeritus_comm_team
    comm_team -= {"reshamas"}  # Exclude reshams from comm_team for web display.

    # Get profiles from GitHub.
    core_devs = [get_profile(login) for login in core_devs]
    emeritus = [get_profile(login) for login in emeritus]
    contributor_experience_team = [get_profile(login) for login in contributor_experience_team]
    emeritus_contributor_experience_team = [get_profile(login) for login in emeritus_contributor_experience_team]
    comm_team = [get_profile(login) for login in comm_team]
    emeritus_comm_team = [get_profile(login) for login in emeritus_comm_team]
    documentation_team = [get_profile(login) for login in documentation_team]

    # Sort teams by last name.
    core_devs = sorted(core_devs, key=key)
    emeritus = sorted(emeritus, key=key)
    contributor_experience_team = sorted(contributor_experience_team, key=key)
    emeritus_contributor_experience_team = sorted(emeritus_contributor_experience_team, key=key)
    documentation_team = sorted(documentation_team, key=key)
    comm_team = sorted(comm_team, key=key)
    emeritus_comm_team = sorted(emeritus_comm_team, key=key)

    return (
        core_devs,
        emeritus,
        contributor_experience_team,
        emeritus_contributor_experience_team,
        comm_team,
        emeritus_comm_team,
        documentation_team,
    )


def get_profile(login):
    """Get the GitHub profile from login."""
    print(f"Getting profile for {login}")
    try:
        profile = get(f"https://api.github.com/users/{login}").json()
    except requests.exceptions.HTTPError:
        return {"name": login, "avatar_url": LOGO_URL, "html_url": ""}

    if profile.get("name") is None:
        profile["name"] = profile.get("login", login)

    # Fix missing names.
    missing_names = {
        "bthirion": "Bertrand Thirion",
        "dubourg": "Vincent Dubourg",
        "Duchesnay": "Edouard Duchesnay",
        "Lars": "Lars Buitinck",
        "MechCoder": "Manoj Kumar",
    }
    if profile["name"] in missing_names:
        profile["name"] = missing_names[profile["name"]]

    return profile


def key(profile):
    """Get a sorting key based on the lower case last name, then first name."""
    components = profile["name"].lower().split()
    return " ".join([components[-1]] + components[:-1])


def generate_table(contributors):
    lines = [
        ".. raw :: html",
        "    <!-- Generated by generate_authors_table.py -->",
        '    <div class="sk-authors-container">',
        "    <style>",
        "      img.avatar { border-radius: 10px; }",
        "    </style>",
    ]
    for contributor in contributors:
        lines.append("    <div>")
        lines.append(
            f"    <a href='{contributor['html_url']}'><img src='{contributor['avatar_url']}' class='avatar' /></a> <br />"
        )
        lines.append(f"    <p>{contributor['name']}</p>")
        lines.append("    </div>")
    lines.append("    </div>")
    return "\n".join(lines) + "\n"


def generate_list(contributors):
    lines = [f"- {contributor['name']}" for contributor in contributors]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    (
        core_devs,
        emeritus,
        contributor_experience_team,
        emeritus_contributor_experience_team,
        comm_team,
        emeritus_comm_team,
        documentation_team,
    ) = get_contributors()

    print("Generating rst files")
    # Write output files using Path.write_text for simplicity.
    (REPO_FOLDER / "doc" / "maintainers.rst").write_text(generate_table(core_devs), encoding="utf-8")
    (REPO_FOLDER / "doc" / "maintainers_emeritus.rst").write_text(generate_list(emeritus), encoding="utf-8")
    (REPO_FOLDER / "doc" / "contributor_experience_team.rst").write_text(
        generate_table(contributor_experience_team), encoding="utf-8"
    )
    (REPO_FOLDER / "doc" / "contributor_experience_team_emeritus.rst").write_text(
        generate_list(emeritus_contributor_experience_team), encoding="utf-8"
    )
    (REPO_FOLDER / "doc" / "communication_team.rst").write_text(generate_table(comm_team), encoding="utf-8")
    (REPO_FOLDER / "doc" / "communication_team_emeritus.rst").write_text(
        generate_list(emeritus_comm_team), encoding="utf-8"
    )
    (REPO_FOLDER / "doc" / "documentation_team.rst").write_text(generate_table(documentation_team), encoding="utf-8")
