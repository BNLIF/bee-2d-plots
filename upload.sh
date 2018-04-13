#!/bin/sh

rsync -avz -e ssh \
--exclude 'Thumbs.db' \
tmp_plots/ \
wirecell@twister.phy.bnl.gov:/data1/lartpc/plots/uboone/
