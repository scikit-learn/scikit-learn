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

The scikit-learn community is organized in **Roles**. The Roles in the project
governance are listed below sorted by increasing involvement in the project:

- Contributors
- Recurring Contributors
- Core Contributors

Roles' Definition, Responsibilities, Rights and Recognition are inherited.

Contributors
------------

Definition
~~~~~~~~~~

**Contributors are individuals or organisations who participate technically, socially
or financially to the life or improvement of the scikit-learn.**

.. note:: **Contributors aren't necessarily developers.**

There are many ways to contribute to the project:

- improving the documentation (with respect to its UX and UI, structural organisation,
  content, examples and references)
- organising events (such as sprints, workshop, meetups, ...)
- reviewing Pull Requests
- producing qualitative reference materials (blog posts, presentation, ...)
- mentoring other Contributors
- triaging Issues
- authoring Pull Requests
- participating in the project's discussions
- advocating the project to the public
- searching for funds for the project
- funding the project

Responsibilities
~~~~~~~~~~~~~~~~

- Must personally uphold the Code of Conduct.

Rights and Recognition
~~~~~~~~~~~~~~~~~~~~~~

- Have their contribution acknowledge somewhere (e.g. listed in the changelog)

Recurring Contributors
----------------------

Definition
~~~~~~~~~~

**Recurring Contributors are Contributors who recurrently participate to the project.**

Recurring Contributors are seen as Contributors who get confortable with the project and
who can provide valuable insight on some of its aspects.

Recurrent Contributor role will often be an intermediate step for people in becoming
Core Contributors once their contributions are frequent enough and during a sustained
period of time if they are willing to get involved further in the project.

Recurrent Contributor is an valuable role by itself for people who want to be part of
the project on a more advisory-like role, as they for example might not have the time
availability or don't want the responsibilities that come with being a Core Contributor.

Contributors are promoted as Recurrent Contributors after the proposal of at least
another Recurrent Contributor and after a vote on the private mailing list.

Responsibilities
~~~~~~~~~~~~~~~~

- Willing to attend events and sessions organised by their team.

Rights and Recognition
~~~~~~~~~~~~~~~~~~~~~~

- Are part of the scikit-learn GitHub organisation
- Are listed as Recurring Contributors on scikit-learn websites
- Can nominate Contributors to be Recurrent Contributors
- Might become Core-Contributors
- Can participate to monthly Core-contributors meetings


Core Contributors
-----------------

Definition
~~~~~~~~~~

**Core Contributors are Recurring Contributors which entrusted with their involvement to
the well being of the Project thanks to their frequency of quality contributions over a
sustained period of time.**

.. note:: **Core Contributors aren't necessarily Core Developers.**


They are the main governing and decision body of the Project and are therefore given
voting and managing rights to the Project services.

Team memberships for Core Contributors refer to their Core Contributor role, as
experienced and entrusted members of the team they are considered Recurrent Contributors
for the rest of the teams and given permissions accordingly.

The exact permissions of all Core Contributors may therefore not be the same and depend
on their team memberships. Even if they have commit rights, Core Contributors should
still have their pull requests reviewed by at least one other Core Contributor
before merging.

In rare circumstances, a Core Contributor may bypass this review if in their best
judgement it is appropriate to do so, but such expedited PR merges must be justifiable
and ultimately subject to review post hoc, if necessary.

Recurrent Contributors are promoted as Core Contributors after the proposal of at least
another Core Contributor and after a vote on the private mailing list.

Responsibilities
~~~~~~~~~~~~~~~~

- Willing to organise events and sessions for their team

Rights and Recognition
~~~~~~~~~~~~~~~~~~~~~~

- Are part of the Maintainers GitHub Team
- Are listed as Core Contributors on scikit-learn websites
- Can nominate Recurrent Contributors to be Core Contributors
- Have “commit rights” to the repository
- Take part in the projects' governance discussions
- Have vote rights for the project direction (promoting a contributor to be a
  core-developer, approving a SLEP, etc.)


Contributors
------------

Contributors are community members who contribute in concrete ways to the
project. Anyone can become a contributor, and contributions can take many forms
– not only code – as detailed in the :ref:`contributors guide <contributing>`.

Contributor Experience Team
---------------------------

The contributor experience team is composed of community members who have permission on
github to label and close issues. :ref:`Their work <bug_triaging>` is
crucial to improve the communication in the project and limit the crowding
of the issue tracker.

Similarly to what has been decided in the `python project
<https://devguide.python.org/triaging/#becoming-a-member-of-the-python-triage-team>`_,
any contributor may become a member of the scikit-learn contributor experience team,
after showing some continuity in participating to scikit-learn
development (with pull requests and reviews).
Any core developer or member of the contributor experience team is welcome to propose a
scikit-learn contributor to join the contributor experience team. Other core developers
are then consulted: while it is expected that most acceptances will be
unanimous, a two-thirds majority is enough.
Every new member of the contributor experience team will be announced in the mailing
list. Members of the team are welcome to participate in `monthly core developer meetings
<https://github.com/scikit-learn/administrative/tree/master/meeting_notes>`_.

.. _communication_team:

Communication team
-------------------

Members of the communication team help with outreach and communication
for scikit-learn. The goal of the team is to develop public awareness of
scikit-learn, of its features and usage, as well as branding.

For this, they can operate the scikit-learn accounts on various social
networks and produce materials.

Every new communicator will be announced in the mailing list.
Communicators are welcome to participate in `monthly core developer meetings
<https://github.com/scikit-learn/administrative/tree/master/meeting_notes>`_.

Core developers
---------------

Core developers are community members who have shown that they are dedicated to
the continued development of the project through ongoing engagement with the
community. They have shown they can be trusted to maintain scikit-learn with
care. Being a core developer allows contributors to more easily carry on
with their project related activities by giving them direct access to the
project's repository and is represented as being an organization member on the
scikit-learn `GitHub organization <https://github.com/orgs/scikit-learn/people>`_.
Core developers are expected to review code
contributions, can merge approved pull requests, can cast votes for and against
merging a pull-request, and can be involved in deciding major changes to the
API.

New core developers can be nominated by any existing core developers. Once they
have been nominated, there will be a vote by the current core developers.
Voting on new core developers is one of the few activities that takes place on
the project's private management list. While it is expected that most votes
will be unanimous, a two-thirds majority of the cast votes is enough. The vote
needs to be open for at least 1 week.

Core developers that have not contributed to the project (commits or GitHub
comments) in the past 12 months will be asked if they want to become emeritus
core developers and recant their commit and voting rights until they become
active again. The list of core developers, active and emeritus (with dates at
which they became active) is public on the scikit-learn website.

Technical Committee
-------------------
The Technical Committee (TC) members are core developers who have additional
responsibilities to ensure the smooth running of the project. TC members are expected to
participate in strategic planning, and approve changes to the governance model.
The purpose of the TC is to ensure a smooth progress from the big-picture
perspective. Indeed changes that impact the full project require a synthetic
analysis and a consensus that is both explicit and informed. In cases that the
core developer community (which includes the TC members) fails to reach such a
consensus in the required time frame, the TC is the entity to resolve the
issue.
Membership of the TC is by nomination by a core developer. A nomination will
result in discussion which cannot take more than a month and then a vote by
the core developers which will stay open for a week. TC membership votes are
subject to a two-third majority of all cast votes as well as a simple majority
approval of all the current TC members. TC members who do not actively engage
with the TC duties are expected to resign.

The Technical Committee of scikit-learn consists of :user:`Thomas Fan
<thomasjpfan>`, :user:`Alexandre Gramfort <agramfort>`, :user:`Olivier Grisel
<ogrisel>`, :user:`Adrin Jalali <adrinjalali>`, :user:`Andreas Müller
<amueller>`, :user:`Joel Nothman <jnothman>`, :user:`Gaël Varoquaux
<GaelVaroquaux>` and :user:`Roman Yurchak <rth>`.

Decision Making Process
=======================
Decisions about the future of the project are made through discussion with all
members of the community. All non-sensitive project management discussion takes
place on the project contributors' `mailing list <mailto:scikit-learn@python.org>`_
and the `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_.
Occasionally, sensitive discussion occurs on a private list.

Scikit-learn uses a "consensus seeking" process for making decisions. The group
tries to find a resolution that has no open objections among core developers.
At any point during the discussion, any core-developer can call for a vote, which will
conclude one month from the call for the vote. Any vote must be backed by a
:ref:`SLEP <slep>`. If no option can gather two thirds of the votes cast, the
decision is escalated to the TC, which in turn will use consensus seeking with
the fallback option of a simple majority vote if no consensus can be found
within a month. This is what we hereafter may refer to as “the decision making
process”.

Decisions (in addition to adding core developers and TC membership as above)
are made according to the following rules:

* **Minor Documentation changes**, such as typo fixes, or addition / correction of a
  sentence, but no change of the scikit-learn.org landing page or the “about”
  page: Requires +1 by a core developer, no -1 by a core developer (lazy
  consensus), happens on the issue or pull request page. Core developers are
  expected to give “reasonable time” to others to give their opinion on the pull
  request if they're not confident others would agree.

* **Code changes and major documentation changes**
  require +1 by two core developers, no -1 by a core developer (lazy
  consensus), happens on the issue of pull-request page.

* **Changes to the API principles and changes to dependencies or supported
  versions** happen via a :ref:`slep` and follows the decision-making process outlined above.

* **Changes to the governance model** use the same decision process outlined above.


If a veto -1 vote is cast on a lazy consensus, the proposer can appeal to the
community and core developers and the change can be approved or rejected using
the decision making procedure outlined above.

.. _slep:

Enhancement proposals (SLEPs)
==============================
For all votes, a proposal must have been made public and discussed before the
vote. Such proposal must be a consolidated document, in the form of a
"Scikit-Learn Enhancement Proposal" (SLEP), rather than a long discussion on an
issue. A SLEP must be submitted as a pull-request to
`enhancement proposals <https://scikit-learn-enhancement-proposals.readthedocs.io>`_
using the `SLEP template <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep_template.html>`_.
