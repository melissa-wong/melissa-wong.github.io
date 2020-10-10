---
title: Projects
permalink: /projects/
---

# MAADCAP
The purpose of this project was to develop a process for prioritizing investment and acquistion decisions for cyber capabilities for a government organization.  The project had previously developed a framework and collected data by the time I joined the team.  The data consisted of subject matter expert (SME) ratings on data sources, threat techniques and analytic capabilities. I tackled the problem of how to quantitatively assess consensus among the SME ratings.

- First, I did a literature search on measures of consensus so that I could recommended [suitable statistics](https://melissa-wong.github.io/doc/MAADCAP-Notes.pdf).

- Next, I looked at how to incoporate self-reported rater confidence data into the analysis; the team had previously collected this information but had no idea how to make use of it. I ended up [adapting a Markov Chain model of the Delphi method](https://melissa-wong.github.io/doc/Using_Confidence.pdf) developed by [DeGroot](http://socsci2.ucsd.edu/~aronatas/project/academic/degroot%20consensus.pdf).

- Finally, once we were able to quantitatively assess consensus, the question was how to priortize data sources and/or analytic investments. For example, the government organization may want to prioritize the minimum set of data sources which give the maximum coverage across threats. This is essentially a [modified set cover problem](https://melissa-wong.github.io/doc/set_cover.pdf).
