# Scalable_SLS

This repository contains code for a scalable System Level Synthesis algorithm for optimal control of output-feedback switched systems. The first idea is to synthesize a controller for each possible k-length subword of the switching signal, given the entire language of possible switching signals. This idea is improved using a clustering approach, where $$N_c$$ clusters are formed based on system Markov parameters $$CA^{i-1}B$$, and one controller is assigned to each cluster.

These efforts seek to improve computational complexity in comparison to our [earlier work which used mode-prefix-based controllers](https://arxiv.org/abs/2505.13105). This repository is under active development, and the code is to be considered incomplete.
