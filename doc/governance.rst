.. _governance:

===========================================
Scikit-learn Governance and Decision-Making
===========================================

The purpose of this document is to formalize the governance process used by the
scikit-learn project, to clarify how decisions are made and how the various
elements of our community interact.
This document establishes a decision-making structure that takes into account
feedback from all members of the community and strives to find consensus, while
avoiding any deadlocks.

This is a meritocratic, consensus-based community project with democratic elements.
Anyone with an interest in the project can join the community, contribute to the
project design and participate in the decision-making process. This document describes
how that participation takes place and how to set about earning merit within
the project community.

Roles And Responsibilities
==========================

There are three basic roles.

Contributors
------------

Contributors are persons who contribute in concrete ways to the project, i.e. they
invest time and energy to improve the project for the greater community.
Anyone can become a contributor, and contributions can take many forms:
coding, reviewing pull request, communications, organization, triaging, replying
to mails in the mailing list, etc.
Every contributor is welcome to participate in `monthly core contributor meetings
<https://github.com/scikit-learn/administrative/tree/master/meeting_notes>`_.
More details can be found in the :ref:`contributors guide <contributing>`.

Core Contributors
-----------------

Core contributors are persons who have shown that they are dedicated to the continued
development of the project through ongoing engagement with the community.
They have shown that they can be trusted to maintain scikit-learn with care.
Core contributors are expected to invest continued time and energy in the project.
This includes selectively (not cumulatively):

  - Ensure the maintenance of scikit-learn
  - Help other contributors to engage in the project and to carry on with their project
    related activities.
  - Conduct public relations work
  - Fundraising
  - Triage github issues and pull requests
  - Review code
  - Contribute code
  - Create a common development roadmap

Therefore, core contributors are given the following rights:

  - Call for votes and cast votes according to the decision-making process as detailed
    below.
  - Be involved in deciding major changes to the API and the governance.
  - Direct access to the project's GitHub repository and joining as organization member
    in the scikit-learn `GitHub organization
    <https://github.com/orgs/scikit-learn/people>`_.
    This right might not be exercised for security reasons, but it may be claimed at
    any time.
  - Cast a vote for, i.e. approve, or cast a vote against pull-requests.
  - Merge pull-requests.

New core contributors can be nominated by any existing core contributor.
Once they have been nominated, there will be a vote by the current core contributors.
Voting on new core contributors is one of the few activities that takes place on the
project's private management list.
While it is expected that most votes will be unanimous, a two-thirds majority of the
cast votes is needed for acceptance of the candidate as new core contributor.
The vote needs to be open for at least 2 and at most 4 weeks.

Core contributors can step back from their role and privileges at any time.
If they have not reasonably contributed to the project in the past 12 months, they will
be kindly asked to step back by the technical committee.
As *ultima ratio*, the technical committee is allowed to call for a vote to withdraw
the core contributor role from a current member.
A 3/4 majority of all current core contributors is needed that this vote passes.
It is open for exactly 4 weeks and is done on a private list.
Core contributors that step back can, if they want to, be listed as emeritus core
contributors on the public website.

The list of core contributors, active and emeritus (with dates at
which they became active) is public on the scikit-learn website.

# TODO: remove/replace .. _communication_team:
# TODO: replace core developer by core contributor in all files.
# TODO: what to do with members of the contributor experience team, communication team?

Technical Committee
-------------------
The Technical Committee (TC) members are core contributors who have additional
responsibilities to ensure the smooth running of the project. TC members are expected to
participate in strategic planning, and approve changes to the governance model.
The purpose of the TC is to ensure a smooth progress from the big-picture
perspective. Indeed changes that impact the full project require a synthetic
analysis and a consensus that is both explicit and informed. In cases that the
core contributor community (which includes the TC members) fails to reach such a
consensus in the required time frame, the TC is the entity to resolve the
issue.
Membership of the TC is by nomination by a core contributor. A nomination will
result in discussion which cannot take more than a month and then a vote by
the core contributors which will stay open for a week. TC membership votes are
subject to a two-third majority of all cast votes as well as a simple majority
approval of all the current TC members. TC members who do not actively engage
with the TC duties are expected to resign.

The Technical Committee of scikit-learn consists of :user:`Thomas Fan
<thomasjpfan>`, :user:`Alexandre Gramfort <agramfort>`, :user:`Olivier Grisel
<ogrisel>`, :user:`Adrin Jalali <adrinjalali>`, :user:`Andreas Müller
<amueller>`, :user:`Joel Nothman <jnothman>` and :user:`Gaël Varoquaux
<GaelVaroquaux>`.

Decision-Making Process
=======================
Decisions about the future of the project are made through discussion with all
members of the community. All non-sensitive project management discussion takes
place on the project contributors' `mailing list <mailto:scikit-learn@python.org>`_
and the `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_.
Occasionally, sensitive discussion occurs on a private list.

Scikit-learn uses a "consensus seeking" process for making decisions. The group
tries to find a resolution that has no open objections among core contributors.
At any point during the discussion, any core contributor can call for a vote, which
will conclude one month from the call for the vote. Any vote must be backed by a
:ref:`SLEP <slep>`. If no option can gather two thirds of the votes cast, the
decision is escalated to the TC, which in turn will use consensus seeking with
the fallback option of a simple majority vote if no consensus can be found
within a month. This is what we hereafter may refer to as "**the decision-making
process**".

Decisions (in addition to adding core contributors and TC membership as above)
are made according to the following rules:

* **Minor Documentation changes**, such as typo fixes, or addition / correction of a
  sentence, but no change of the scikit-learn.org landing page or the “about”
  page: Requires +1 by a core contributor, no -1 by a core contributor (lazy
  consensus), happens on the issue or pull request page. Core contributors are
  expected to give “reasonable time” to others to give their opinion on the pull
  request if they're not confident others would agree.

* **Code changes and major documentation changes**
  require +1 by two core contributors, no -1 by a core contributor (lazy
  consensus), happens on the issue of pull-request page.

* **Changes to the API principles and changes to dependencies or supported
  versions** happen via a :ref:`slep` and follows the decision-making process outlined above.

If a veto -1 vote is cast on a lazy consensus, the proposer can appeal to the
community and core contributors and the change can be approved or rejected using
the decision-making procedure outlined above.

Governance Model Changes
------------------------

Governance model changes occur through an enhancement proposal or a GitHub Pull
Request. An enhancement proposal will go through "**the decision-making process**"
described in the previous section. Alternatively, an author may propose a change
directly to the governance model with a GitHub Pull Request. Logistically, an
author can open a Draft Pull Request for feedback and follow up with a new
revised Pull Request for voting. Once that author is happy with the state of the
Pull Request, they can call for a vote on the public mailing list. During the
one-month voting period, the Pull Request can not change. A Pull Request
Approval will count as a positive vote, and a "Request Changes" review will
count as a negative vote. If two-thirds of the cast votes are positive, then
the governance model change is accepted.

.. _slep:

Enhancement proposals (SLEPs)
==============================
For all votes, a proposal must have been made public and discussed before the
vote. Such proposal must be a consolidated document, in the form of a
"Scikit-Learn Enhancement Proposal" (SLEP), rather than a long discussion on an
issue. A SLEP must be submitted as a pull-request to
`enhancement proposals <https://scikit-learn-enhancement-proposals.readthedocs.io>`_
using the `SLEP template <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep_template.html>`_.
