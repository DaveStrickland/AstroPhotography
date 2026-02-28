.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are appreciated! Every little bit helps, and 
credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Contributions need not require any coding skill. Bug reports or suggestions
for enhancements are valuable in themselves. If you do wish to contribute
to code or documentation please familiarize yourself with the [AI Policy](AI_POLICY.md).
This project places limits on the forms of AI assistance allowed to reduce
the burden AI slop places on project maintainers.

Report Bugs or Suggest New Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Report bugs or suggest new features or enhancements at 
https://github.com/DaveStrickland/astrophotography/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Explain what you were trying to do, what files you were using, what happened, and what
  you expected to happen.
* Providing the steps needed to reproduce the problem greatly increases the chance
  of the bug report being worked and fixed. One-off weird behavior that can't be
  reproduced will likely never be addressed. 
* Log files, in particular running the code with `--loglevel=DEBUG`, are invaluable
  in diagnosting bugs, so please include or attach them to the Issues.
* If the data file (or a few files) are not too large, consider attaching them or 
  providing a Dropbox link where they can be downloaded.

You can also suggest new features, or enhancements to existing features.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Any issue tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

AstroPhotography could always use more documentation, whether as part of the
official AstroPhotography docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/DaveStrickland/astrophotography/issues.

If you are proposing a new feature:

* Explain in detail what the need is, and how it might be implemented.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, contributions are welcome, but the
  pace of development (and response to new issues, PRs, comments, etc) is
  slow. 

Get Started With Local Development!
-----------------------------------

Ready to contribute by coding? Here's how to set up `astrophotography` for local development.

1. Fork the `astrophotography` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/astrophotography.git

3. Follow the basic setup and installation instructions in the [README](README.md).

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   It is recommended that branches are named after the issue being worked. If you are working
   on Issue 66 you should name the branch `issue-066`. Now you can make your changes locally.

5. When you're done making changes, check that your changes pass any tests you may have
   created youself, e.g. by running a script or set of calls against some input data
   files.

   Follow the Developers Only and Documentation instructions in the [README](README.md).
   The test suite is extremely limited at present. The documentation is in a better
   state, and you should check that the documentation for files and functions you have
   altered is up-to-date and correct.

6. Commit your changes and push your branch to GitHub::

    $ git add xxx # altered or added files
    $ git commit -m "issue-xxx Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

    Any required AI disclosure to the commit should be done at this time. See the 
    [AI Policy](AI_POLICY.md).

7. Submit a pull request through the GitHub website. The Pull Request template will lead you
   through the information you will need to provide.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines. You
can see what is required by looking at the [pull_request_template.md](.github/PULL_REQUEST_TEMPLATE/pull_request_template.md)

Pull Request titles should reflect the name of the Issue they are associated with.
Pull Requests should *not* attempt to automatically close an Issue using text of the form
`Closes #66` or `Fixes #66`. It is up to the maintainers to assess whether a contribution
fully addresses the Issue or whether other changes are necessary.

Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in RELEASE_NOTES.rst).
Then run::

$ bump2version patch # possible: major / minor / patch
$ git push
$ git push --tags

