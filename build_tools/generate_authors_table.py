"""
This script generates an html table of contributors, with names and avatars.
The list is generated from scikit-learn's teams on GitHub, plus a small number
of hard-coded contributors.

The table should be updated for each new inclusion in the teams.
Generating the table requires admin rights.
"""
import sys
import requests
import getpass
import time

print("user:", file=sys.stderr)
user = input()
passwd = getpass.getpass("Password or access token:\n")
auth = (user, passwd)

LOGO_URL = 'https://avatars2.githubusercontent.com/u/365630?v=4'

def get(url):
    for sleep_time in [10, 30, 0]:
        reply = requests.get(url, auth=auth)
        api_limit = ("message" in reply.json()
                     and "API rate limit exceeded" in reply.json()["message"])
        if not api_limit:
            break
        print("API rate limit exceeded, waiting..")
        time.sleep(sleep_time)

    reply.raise_for_status()
    return reply


def group_iterable(iterable):
    """Group iterable into lines"""
    group = []
    for element in iterable:
        group.append(element)
    if len(group) != 0:
        yield group


def get_contributors():
    """Get the list of contributor profiles. Require admin rights."""
    # get members of scikit-learn core-dev on GitHub
    core_devs = []
    team = 11523
    for page in [1, 2]:  # 30 per page
        reply = get("https://api.github.com/teams/%d/members?page=%d" %
                    (team, page))
        core_devs.extend(reply.json())

    # get members of scikit-learn on GitHub
    members = []
    for page in [1, 2]:  # 30 per page
        reply = get(
            "https://api.github.com/orgs/scikit-learn/members?page=%d" %
            (page, ))
        members.extend(reply.json())

    # keep only the logins
    core_devs = [c['login'] for c in core_devs]
    members = [c['login'] for c in members]

    # add missing contributors with GitHub accounts
    members.extend(['dubourg', 'mbrucher', 'thouis', 'jarrodmillman'])
    # add missing contributors without GitHub accounts
    members.extend(['Angel Soler Gollonet'])
    # remove CI bots
    members.remove('sklearn-ci')
    members.remove('sklearn-lgtm')
    members.remove('sklearn-wheels')

    # remove duplicate, and get the difference of the two sets
    core_devs = set(core_devs)
    members = set(members)
    emeritus = members.difference(core_devs)

    # get profiles from GitHub
    core_devs = [get_profile(login) for login in core_devs]
    emeritus = [get_profile(login) for login in emeritus]

    # sort by last name
    core_devs = sorted(core_devs, key=key)
    emeritus = sorted(emeritus, key=key)

    return core_devs, emeritus


def get_profile(login):
    """Get the GitHub profile from login"""
    print("get profile for %s" % (login, ))
    try:
        profile = get("https://api.github.com/users/%s" % login).json()
    except requests.exceptions.HTTPError:
        return dict(name=login, avatar_url=LOGO_URL, html_url="")

    if profile["name"] is None:
        profile["name"] = profile["login"]

    # fix missing names
    missing_names = {
        'bthirion': 'Bertrand Thirion',
        'dubourg': 'Vincent Dubourg',
        'Duchesnay': 'Edouard Duchesnay',
        'Lars': 'Lars Buitinck',
        'MechCoder': 'Manoj Kumar',
        'jeremiedbb': 'Jérémie Du Boisberranger',
    }
    if profile["name"] in missing_names:
        profile["name"] = missing_names[profile["name"]]

    return profile


def key(profile):
    """Get the last name in lower case"""
    return profile["name"].split(' ')[-1].lower()


def generate_table(contributors):
    lines = [
        (".. raw :: html\n"),
        ("    <!-- Generated by generate_authors_table.py -->"),
        ("    <div class='sk-authors-container'>"),
        ("    <style>"),
        ("      img.avatar {border-radius: 10px;}"),
        ("    </style>"),
    ]
    for contributor in group_iterable(contributors):
        lines.append("    <div>")
        lines.append(
            "    <a href='%s'><img src='%s' class='avatar' /></a> <br />" %
            (contributor["html_url"], contributor["avatar_url"]))
        lines.append("    <p>%s</p>" % (contributor["name"], ))
        lines.append("    </div>")
    lines.append("    </div>")
    return '\n'.join(lines)


def generate_list(contributors):
    lines = []
    for contributor in contributors:
        lines.append("- %s" % (contributor["name"], ))
    return '\n'.join(lines)


if __name__ == "__main__":

    core_devs, emeritus = get_contributors()

    with open("../doc/authors.rst", "w+") as rst_file:
        rst_file.write(generate_table(core_devs))

    with open("../doc/authors_emeritus.rst", "w+") as rst_file:
        rst_file.write(generate_list(emeritus))
