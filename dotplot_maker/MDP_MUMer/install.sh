#!/usr/bin/env bash

make
if ! [[ -d ../bin ]]; then
    mkdir ../bin
fi
cd ../bin
ln -s -T ../MDP_MUMer/mummer MDP_mummer
ln -s -T ../MDP_MUMer/nucmer MDP_nucmer
ln -s -T ../MDP_MUMer/mummerplot MDP_mummerplot
ln -s -T ../MDP_MUMer/delta-filter MDP_delta-filter
ln -s -T ../MDP_MUMer/show-coords MDP_show-coords

