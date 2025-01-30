How is the GEOS project organized
=================================

The **Architecture Committee** oversees significant changes to the GEOS architecture, ensuring compliance with organizational standards and security measures. Developers must present their ideas to the committee for initial approval before opening a Pull Request (PR), with the committee aiming to provide decisions within a week. The committee meets bi-weekly and on-demand, including representatives from each entity and relevant stakeholders.

The **Operations Committee** handles the day-to-day management of the project, including managing pull requests (PRs), resolving tickets/issues, and assigning tasks. This committee ensures smooth development workflows and aligns daily operations with project goals.

Tools:
- Everything currently happens on GitHub.
- New tools will be evaluated depending on the needs.


Architecture Committee
----------------------

## 1. General information

- **Committee Overview:**
  - The architecture committee overseea changes to the system architecture.
  - Prior to opening a Pull Request (PR), developers must present their idea to the committee for initial feedback and approval.

- **Meeting Schedule:**
  - The committee meets bi-weekly and on-demand.
  - Developers must request a slot for their presentation via the ``#dev-infrastructure`` Slack channel.
    If you do not have access to our Slack workspace, you can open an issue on GitHub requesting a review.

- **Committee Members:**
  - The committee consists of one representative for partner organization.
  - Additional members with relevant skills and the product owner of affected components are included.
  - Meetings are open for all interested parties to attend.

- **Changes Requiring Committee Review:**
  - Any development significantly modifying the architecture.
  - Changes involving logs, ghosting, adding/removing TPL dependencies, or high-order methods (e.g., Shiva).
  - Significant modifications to the architecture must be validated by the committee before merging.

- **Changes Not Requiring Committee Review:**
  - Bug fixes, documentation updates, and numerical methods tuning.

## 2. Presenting your idea to the committee
- **Description of Change:**
  - Describe the proposed change, including the motivation and problem it solves.
  - Consider making this proposal using markdown or plantuml to include in the documentation when merged.

- **Impact Analysis:**
  - Analyze and document the impact on the existing system, including dependencies, affected modules, and performance considerations.

- **Implementation Details:**
  - Provide technical details of the implementation, including diagrams if necessary.

- **Testing Strategy:**
  - Outline a comprehensive testing plan, including unit tests, integration tests, and other relevant testing methods.

- **Initial Presentation:**
  - The committee aims to provide a go/no-go decision within one week.

## 3. Implementation and final approval
- **Approval Notification:**
  - Await official approval notification from the architecture committee.
- **Merge Changes:**
  - Once approved, the standard PR process applies.


Operations Committee
--------------------

## 1. General information
- **Primary Responsibilities:**
  - Oversee the review and merging of pull requests.
  - Manage and triage tickets/issues.
  - Assign tasks to developers and track their progress.

- **Meeting Schedule:**
  - The committee meets weekly to review progress, manage the current backlog, and plan upcoming work. Frequency may decrease as the backlog diminishes.
  - Additional meetings may be scheduled as needed for urgent issues.

- **Committee Members:**
  - Includes representatives from key areas of the project, such as team leads, senior developers, and project managers.
  - Meetings are open to other team members involved in specific tasks or issues being discussed.

## 2. Pull Request (PR) Management
- **Reviewer Selection:**
  - Confirm the selection of reviewers for all developments approved by the **Architecture Committee**.
  - Ensure the _module owners_ of the impacted components of GEOS are included in the review process.
  - Involve a mix of junior and senior developers in each PR to promote knowledge sharing and thorough reviews.

- **Review Process:**
  - Reviewers ensure that code meets coding standards, passes all tests, and does not introduce new issues.
  - Provide feedback through PR comments, and developers address any concerns before the PR can be merged.

- **Merging PRs:**
  - Once a PR is approved by the required number of reviewers, it can be merged into the main branch.
  - While one validation is required, multiple approvals are recommended for major changes.

## 3. Issue Management
- **Triaging:**
  - The committee triages incoming issues, prioritizing on urgency and impact.
  - The committee updates the status of tickets to ensure high-priority issues are addressed promptly.

- **Meeting Frequency:**
  - Meet weekly to manage the current backlog and discuss issue resolution. Frequency may be adjusted based on backlog size.
  - Track the progress of each issue to ensure timely resolution.

- **Progress Monitoring:**
  - Regularly check on the progress of assigned tasks.
  - Address any blockers or issues that may prevent tasks from being completed on time.

## Notes
- Continuous improvement and adaptation are key to addressing evolving project needs and challenges.
- For parties involved in the _FC Maelstrom_ project:
  - Each party can assign tasks to their team members.
  - In cases where resources need to move between parties, collegial discussions will occur during the weekly meetings.
