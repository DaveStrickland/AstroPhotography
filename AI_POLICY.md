# AI-Assisted Contributions Policy

This is the AI Policy **for** _human contributors_ about **using** AI 
that describes:

- What is allowed and what is disallowed
- Disclosure obligations

For instructions **to** _AI assistants_ about **providing** AI assistance, 
please see [AGENTS](AGENTS.md).

In the following document the term AI is used to refer to
generative-AI tools (LLMs) of the form popularized by ChatGPT. 
It does not refer to machine learning (ML) in general.
AI-tools are tools such as Microsoft Copilot, Cursor, and others
using LLM-based models such as ChatGPT, Claude, and others.
Non-AI tools are such as IntelliSense, PyLance, ruff, and the like.

## Summary

- If possible, avoid the use of AI altogether.
- AI tool use that modifies any part of this project should be 
  limited, targeted, and human directed. 
- You must read, review, and understand every AI-based suggestion
  or modification, and be prepared to explain what it does and why
  it is necessary to another human.
- AI shall not be used to modify core logic AND write or update the
  test for that logic.
- Adding or updating docstrings is a valid low-risk use of AI tools.
- Bots, including contributions from bot accounts or accounts that 
  look like bots, are prohibited.
- AI shall not be used to create or author Issues or Pull Requests
  directly. 
- You must disclose the use of AI tools that have resulted in modified 
  files, beyond cases equivalent to existing non-AI code assist
  tools (e.g. Intellisense).
- AI tool usage that does not modify files, e.g. "Explain what the
  selected code does" or "Analyze the repository and produce a list
  of files and functions lacking docstrings" cannot be prevented,
  and so it is technically permissable. But remember there are other ways
  to accomplish such tasks that do not involve such ethically dubious,
  environmentally unfriendly, and error prone tools.
- All contributions must meet the project’s standards for inclusion. 

These rules imply Chat-like work flows. Agentic work flows are not
permitted. 

The guidelines given here apply even if you subsequently edit, alter,
or extend an AI-generated changes.

These rules apply to all people without disabilities who have
some fluency in English. People with disabilities, and/or people with
limited fluency in English **who require the aid of an AI** may
use it to provide a larger fraction of their contribution. This would
be assessed on a case by case basis. 

Submissions that violate these rules, or could plausibly be construed
to violate them, will be rejected without comment or recourse.

*Large scale initiatives:* The policy does not cover possible large 
scale initiatives which may significantly change the ways the project 
is structured. Such initiatives need to be discussed separately with 
the Project Lead.

## Mandatory Disclosure and Human Accountability

### Accountability

You MUST take the responsibility for your contribution. Contributing 
to this project means vouching for the quality, license compliance,
provenance, and utility of your submission. You, as a human contributor
are always the author and are fully accountable for the entirety of 
these contributions.

When you submit an Issue or submit code in a PR you are signing your 
name to it and warranting that you, a human, are the **Author**. You 
are also responsible for it, it reflects your judgement and skill. 
The project Pull Request templates will also have AI use attestation field.

As already stated, you must read and understand every part of the changes
you make, whether hand-written or AI-generated. If you don't understand
it yourself, don't submit it.

An AI can never be credited as an author or co-author (it cannot 
hold or transfer copyright). This is one of the reasons that AI-authored
and/or AI-submitted PRs are not acceptable.

### Transparency and Disclosure

AI usage above levels described below must be disclosed. 

You MUST disclose the use of AI tools when AI is used to make significant
changes, as described below, even if you subsequently edit or alter
the AI-generated changes. You SHOULD disclose the other uses of non-coding
AI tools, where it might be useful. Routine use of assistive tools for 
correcting grammar and spelling, or for clarifying language, does not 
require disclosure.

Information about the use of AI tools will help us evaluate their impact, 
build new best practices and adjust existing processes.

For contributions tracked in git, the recommended method is an Assisted-by: 
commit message trailer. For other contributions (e.g. Issue or PR text), 
disclosure should be made in the relevant contribution itself.

Examples:

    Assisted-by: Microsoft Copilot 365

    Assisted-by: Microsoft Copilot using Claude Haiku 4.5

It is preferable to provide both the type and version of the AI model being used,
but in the absence of established industry disclosure standards,
just try your best.

On the command line this can be done with a second `-m` flag on the `git commit`, e.g.
```bash
git commit -m "issue-2001 Investigate monolith" -m "Assisted-by: HAL9000"
```
Commits with the `Assisted-by` tag can be found using `git log --all --grep="Assisted-by" --pretty=format:"%h %s"`

The `Assisted-by` lines should be preserved and rolled up together at any
git `squash`.

Failure to disclose AI use can result in summary PR rejection and/or merged
changes being reverted. 

Contribution & Community *Evaluation*: AI tools may be used to assist human 
*reviewers* by providing analysis and suggestions. You MUST NOT use AI as the 
sole or final arbiter in making a substantive or subjective judgment on 
a contribution. Evaluative AI use should also be disclosed.

#### When in the process should AI use be disclosed?

At latest, when a PR is submitted the git commits should have any necessary `Assisted-by`
messages added. It is preferrable to get into the habit of using `Assisted-by`
when you make commits.

## AI Use Cases

The following list is intended to provide some examples of what you can and cannot do,
and whether doing so requires disclosure. It is not intended as an exhaustive list,
but to be illustrative. There is some degree of subjectivity in these guidelines

### Allowed To Some Degree

- Code assist: Anything of the same form and scope as non-AI code assist tools,
  including spell checking, variable name completion, brace and parenthesis addition
  and checking, syntax error highlighting, or find and replace. ALLOWED, NO DISCLOSURE NEEDED.
  (Unfortunately this is needed because the commercial AI-vendors have strong incentives 
  to disable existing tools to force the adoption and use of their AI tools.) 
- Inline suggestions:
  - A one-line suggestion, for example completing a function call with all 
    contextually-correct variable names or adding a single variable to a docstring. 
    ALLOWED, NO DISCLOSURE. 
  - A multi-line suggestion, for example multiple function calls, for loop or
    if else block, a new docstring block, or many lines of altered docstring. 
    ALLOWED, DISCLOSURE REQUIRED even if subsequently hand-editted.
- AI-review of code or code changes that did not change files: ALLOWED, NO DISCLOSURE. 
- Any AI generated multi-line change that does not violate any of the existing rules.
  ALLOWED, DISCLOSURE REQUIRED
- Boilerplate for a new core class or command line script: The skeleton of a class 
  or function, but does *not* implement the core logic. ALLOWED, DISCLOSURE REQUIRED.
- Boilerplate for a test class: The skeleton of a test class that tests all the 
  functionality of a core implementation class, but does *not* implement the test
  logic. ALLOWED, DISCLOSURE REQUIRED.
- Refactorings by an AI-tool:
  - Refactorings that affect 10 lines or less: ALLOWED, DISCLOSURE REQUIRED.
  - Refactorings that change more than 10 lines but less than 50 lines: Discouraged
    but permissible if you, the author, very carefully review and test it before submitting 
    the change. ALLOWED, DISCLOSURE REQUIRED.
  - Refactorings that change > 50 lines: DISALLOWED.
- AI-generated logic for a single new function: This is where good judgement
  and experience comes into play, so there is no hard and fast absolute rule.
  The smaller the function and/or the simpler the task being accomplished the better, but
  as a general rule it is OK if it is less than 10 logical lines of code, disallowed
  if it is so large it doesn't completely fit in an editor window at one time (without scrolling).
  At an intermediate size of 10-50 lines it is allowed but discouraged. Again you, the author, must
  carefully review and test it before submitting the change. So for changes less than 50 lines it
  is ALLOWED, DISCLOSURE REQUIRED.
- AI-generated changes that affect multiple functions or files at a single time.
  - If only altering docstrings, these are ALLOWED, DISCLOSURE REQUIRED.
  - If only altering type-hinting, these are ALLOWED, DISCLOSURE REQUIRED.

### Banned

- Files created, renamed, or deleted by an AI tool: DISALLOWED.
- Large scale refactoring by an AI tool: Any AI-tool refactoring of more than 50 lines. DISALLOWED.
- New external dependencies added *or suggested* by an AI tool: DISALLOWED.
- New class or function over ~50 lines implemented by an AI tool: DISALLOWED.

Even if allowed, changes may be rejected when under review, just as a purely human-written
change might be rejected or require rework.

## Testing Requirements

Currently this project lacks a testing framework that would allow us to impose
testing requirements on any contributed code, let alone AI-assissted contributions.
At a later date this deficiency will be addressed and testing requirements imposed. 

## Inspiration and Credit

This policy is derived from the `autobahn-python` [AI policy](https://github.com/wamp-proto/wamp-ai/blob/bfb4804ae8fda35db61d2a98821a781b2969c59d/AI_POLICY.md) and the [Fedora AI policy](https://docs.fedoraproject.org/en-US/council/policy/ai-contribution-policy/).