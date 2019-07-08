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

# With authentication: up to 5000 requests per hour.
print("user:", file=sys.stderr)
user = input()
passwd = getpass.getpass("Password or access token:\n")
auth = (user, passwd)

ROW_SIZE = 7
LOGO_URL = 'https://avatars2.githubusercontent.com/u/365630?v=4'


def group_iterable(iterable, size):
    """Group iterable into lines"""
    group = []
    for element in iterable:
        group.append(element)
        if len(group) == size:
            yield group
            group = []
    if len(group) != 0:
        yield group


def get_contributors():
    """Get the list of contributor profiles. Require admin rights."""
    # get members of scikit-learn teams on GitHub
    members = []
    team = 11523
    for page in [1, 2]:  # 30 per page
        reply = requests.get(
            "https://api.github.com/teams/%d/members?page=%d"
            % (team, page), auth=auth)
        reply.raise_for_status()
        members.extend(reply.json())

    # keep only the logins
    logins = [c['login'] for c in members]
    # remove duplicate
    logins = set(logins)

    # get profiles from GitHub
    profiles = [get_profile(login) for login in logins]
    # sort by last name
    profiles = sorted(profiles, key=key)

    return profiles


def get_profile(login):
    """Get the GitHub profile from login"""
    profile = requests.get("https://api.github.com/users/%s" % login,
                           auth=auth).json()
    if 'name' not in profile:
        # default profile if the login does not exist
        return dict(name=login, avatar_url=LOGO_URL, html_url="")
    else:
        if profile["name"] is None:
            profile["name"] = profile["login"]

        # fix missing names
        missing_names = {'bthirion': 'Bertrand Thirion',
                         'Duchesnay': 'Edouard Duchesnay',
                         'Lars': 'Lars Buitinck',
                         'MechCoder': 'Manoj Kumar'}
        if profile["name"] in missing_names:
            profile["name"] = missing_names[profile["name"]]
        return profile


def key(profile):
    """Get the last name in lower case"""
    return profile["name"].split(' ')[-1].lower()


contributors = get_contributors()

print(".. raw :: html\n")
print("    <!-- Generated by generate_authors_table.py -->")
print("    <table>")
print("    <col style='width:%d%%' span='%d'>"
      % (int(100 / ROW_SIZE), ROW_SIZE))
print("    <style>")
print("      img.avatar {border-radius: 10px;}")
print("      td {vertical-align: top;}")
print("    </style>")
for row in group_iterable(contributors, size=ROW_SIZE):
    print("    <tr>")
    for contributor in row:
        print("    <td>")
        print("    <a href='%s'><img src='%s' class='avatar' /></a> <br />"
              % (contributor["html_url"], contributor["avatar_url"]))
        print("    <p>%s</p>" % contributor["name"])
        print("    </td>")
    print("    </tr>")
print("    </table>")
