---
title: Projects
permalink: /projects/
toc: true
---

# MAADCAP
The purpose of this project was to develop a process for prioritizing investment and acquistion decisions for cyber capabilities for a government organization.  The project had previously developed a framework and collected data by the time I joined the team.  The data consisted of subject matter expert (SME) ratings on data sources, threat techniques and analytic capabilities. I tackled the problem of how to quantitatively assess consensus among the SME ratings.

- First, I did a literature search on measures of consensus so that I could recommended [suitable statistics](https://melissa-wong.github.io/doc/MAADCAP-Notes.pdf){:target="_blank" rel="noopener"}.

- Next, I looked at how to incoporate self-reported rater confidence data into the analysis; the team had previously collected this information but had no idea how to make use of it. I ended up [adapting a Markov Chain model of the Delphi method](https://melissa-wong.github.io/doc/Using_Confidence.pdf){:target="_blank" rel="noopener"} developed by [DeGroot](http://socsci2.ucsd.edu/~aronatas/project/academic/degroot%20consensus.pdf).

- Finally, once we were able to quantitatively assess consensus, the question was how to priortize data sources and/or analytic investments. For example, the government organization may want to prioritize the minimum set of data sources which give the maximum coverage across threats. This is essentially a [modified set cover problem](https://melissa-wong.github.io/doc/set_cover.pdf){:target="_blank" rel="noopener"}.

# Data Privacy Consulting Project

As part of my M.S. degree, I worked on a consulting project for a real estate software company. The company was considering making summary data available to their customers, but they wanted to ensure that by doing so they wouldn't inadvertently compromise the private data of any individual brokerage. We investigated several differential privacy mechanisms and made recommendations on how the company could implement an algorithm that would provide a measurable risk of privacy loss. The details are in this [report](https://melissa-wong.github.io/doc/STAT684_Report.pdf){:target="_blank" rel="noopener"}.
