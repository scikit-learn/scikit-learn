.. _governance:

===========================================
Scikit-learn governance and decision-making
===========================================

The purpose of this document is to formalize the governance process used by the
scikit-learn project, to clarify how decisions are made and how the various
elements of our community interact.
This document establishes a decision-making structure that takes into account
feedback from all members of the community and strives to find consensus, while
avoiding any deadlocks.

This is a meritocratic, consensus-based community project. Anyone with an
interest in the project can join the community, contribute to the project
design and participate in the decision making process. This document describes
how that participation takes place and how to set about earning merit within
the project community.

Roles And Responsibilities
==========================

We distinguish between contributors, core contributors, and the technical
committee. A key distinction between them is their voting rights: contributors
have no voting rights, whereas the other two groups all have voting rights,
as well as permissions to the tools relevant to their roles.

Contributors
------------

Contributors are community members who contribute in concrete ways to the
project. Anyone can become a contributor, and contributions can take many forms
– not only code – as detailed in the :ref:`contributors guide <contributing>`.
There is no process to become a contributor: once somebody contributes to the
project in any way, they are a contributor.

Core Contributors
-----------------

All core contributor members have the same voting rights and right to propose
new members to any of the roles listed below. Their membership is represented
as being an organization member on the scikit-learn `GitHub organization
<https://github.com/orgs/scikit-learn/people>`_.

They are also welcome to join our `monthly core contributor meetings
<https://github.com/scikit-learn/administrative/tree/master/meeting_notes>`_.

New members can be nominated by any existing member. Once they have been
nominated, there will be a vote by the current core contributors. Voting on new
members is one of the few activities that takes place on the project's private
mailing list. While it is expected that most votes will be unanimous, a
two-thirds majority of the cast votes is enough. The vote needs to be open for
at least 1 week.

Core contributors that have not contributed to the project, corresponding to
their role, in the past 12 months will be asked if they want to become emeritus
members and recant their rights until they become active again. The list of
members, active and emeritus (with dates at which they became active) is public
on the scikit-learn website.

The following teams form the core contributors group.


Contributor Experience Team
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The contributor experience team improves the experience of contributors by
helping with the triage of issues and pull requests, as well as noticing any
repeating patterns where people might struggle, and to help with improving
those aspects of the project.

To this end, they have the required permissions on github to label and close
issues. :ref:`Their work <bug_triaging>` is crucial to improve the
communication in the project and limit the crowding of the issue tracker.

.. _communication_team:

Communication team
~~~~~~~~~~~~~~~~~~

Members of the communication team help with outreach and communication
for scikit-learn. The goal of the team is to develop public awareness of
scikit-learn, of its features and usage, as well as branding.

For this, they can operate the scikit-learn accounts on various social networks
and produce materials. They also have the required rights to our blog
repository and other relevant accounts and platforms.

Documentation team
~~~~~~~~~~~~~~~~~~

Members of the documentation team engage with the documentation of the project
among other things. They might also be involved in other aspects of the
project, but their reviews on documentation contributions are considered
authoritative, and can merge such contributions.

To this end, they have permissions to merge pull requests in scikit-learn's
repository.

Maintainers
~~~~~~~~~~~

Maintainers are community members who have shown that they are dedicated to the
continued development of the project through ongoing engagement with the
community. They have shown they can be trusted to maintain scikit-learn with
care. Being a maintainer allows contributors to more easily carry on with their
project related activities by giving them direct access to the project's
repository. Maintainers are expected to review code contributions, merge
approved pull requests, cast votes for and against merging a pull-request,
and to be involved in deciding major changes to the API.

Technical Committee
-------------------

The Technical Committee (TC) members are maintainers who have additional
responsibilities to ensure the smooth running of the project. TC members are
expected to participate in strategic planning, and approve changes to the
governance model. The purpose of the TC is to ensure a smooth progress from the
big-picture perspective. Indeed changes that impact the full project require a
synthetic analysis and a consensus that is both explicit and informed. In cases
that the core contributor community (which includes the TC members) fails to
reach such a consensus in the required time frame, the TC is the entity to
resolve the issue. Membership of the TC is by nomination by a core contributor.
A nomination will result in discussion which cannot take more than a month and
then a vote by the core contributors which will stay open for a week. TC
membership votes are subject to a two-third majority of all cast votes as well
as a simple majority approval of all the current TC members. TC members who do
not actively engage with the TC duties are expected to resign.

The Technical Committee of scikit-learn consists of :user:`Thomas Fan
<thomasjpfan>`, :user:`Alexandre Gramfort <agramfort>`, :user:`Olivier Grisel
<ogrisel>`, :user:`Adrin Jalali <adrinjalali>`, :user:`Andreas Müller
<amueller>`, :user:`Joel Nothman <jnothman>` and :user:`Gaël Varoquaux
<GaelVaroquaux>`.

Decision Making Process
=======================
Decisions about the future of the project are made through discussion with all
members of the community. All non-sensitive project management discussion takes
place on the project contributors' `mailing list <mailto:scikit-learn@python.org>`_
and the `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_.
Occasionally, sensitive discussion occurs on a private list.

Scikit-learn uses a "consensus seeking" process for making decisions. The group
tries to find a resolution that has no open objections among core contributors.
At any point during the discussion, any core contributor can call for a vote,
which will conclude one month from the call for the vote. Most votes have to be
backed by a :ref:`SLEP <slep>`. If no option can gather two thirds of the votes
cast, the decision is escalated to the TC, which in turn will use consensus
seeking with the fallback option of a simple majority vote if no consensus can
be found within a month. This is what we hereafter may refer to as "**the
decision making process**".

Decisions (in addition to adding core contributors and TC membership as above)
are made according to the following rules:

* **Minor Documentation changes**, such as typo fixes, or addition / correction
  of a sentence, but no change of the ``scikit-learn.org`` landing page or the
  “about” page: Requires +1 by a maintainer, no -1 by a maintainer (lazy
  consensus), happens on the issue or pull request page. Maintainers are
  expected to give “reasonable time” to others to give their opinion on the
  pull request if they're not confident others would agree.

* **Code changes and major documentation changes**
  require +1 by two maintainers, no -1 by a maintainer (lazy
  consensus), happens on the issue of pull-request page.

* **Changes to the API principles and changes to dependencies or supported
  versions** happen via a :ref:`slep` and follows the decision-making process
  outlined above.

* **Changes to the governance model** follow the process outlined in `SLEP020
  <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep020/proposal.html>`__.

If a veto -1 vote is cast on a lazy consensus, the proposer can appeal to the
community and maintainers and the change can be approved or rejected using
the decision making procedure outlined above.

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
issue. A SLEP must be submitted as a pull-request to `enhancement proposals
<https://scikit-learn-enhancement-proposals.readthedocs.io>`_ using the `SLEP
template
<https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep_template.html>`_.
`SLEP000
<https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep000/proposal.html>`__
describes the process in more detail.
