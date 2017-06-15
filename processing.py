#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Z.-L. Deng
# @Email: dawnmsg@gmail.com
# @Date:   2015-10-23 18:21:59
# @Last Modified by:   zde14
# @Last Modified time: 2015-11-13 16:59:40

import pandas as pd

class Data:

    def __init__(self, data_file, dataset='training', subclass=False):
        self.dataset = dataset
        self.data_file = data_file
        self.df = pd.read_csv(
            self.data_file,
            delim_whitespace=True
        )
        self.df['dataset'] = dataset
        if subclass:
            self.df.set_index(['samples', 'dataset', 'class', 'subclass'], inplace=True)
        else:
            self.df.set_index(['samples', 'dataset', 'class'], inplace=True)

# percentage
def percentage(df):
    data = df.T
    data_norm = 100 * data / data.sum()
    data_norm = data_norm.T

    return(data_norm)


# z-score
def standardization(df, method='max_min'):
    data = df
    if method == 'max_min':
        data_st = (data - data.min()) / (data.max() - data.min())

    elif method == 'z_score':
        data_st = (data - data.mean()) / data.std()

    return(data_st)

