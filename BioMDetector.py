#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Z.-L. Deng
# @Email: dawnmsg@gmail.com
# @Date:   2015-10-23 15:31:54
# @Last Modified by:   zde14
# @Last Modified time: 2015-11-16 15:44:43

import argparse
import processing as dp
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from numpy import arange

# biomarker gene detector based on abundance table
class FS:
    """
    Usage::

    >>> fs = FS(training_stand, training_label)
    """

    def __init__(self, data, label):
        self.data = data
        self.label = label
        self.genes = list(self.data.columns)

    # differential expression gene detection
    def edgeR(self, fdr_cutoff=0.05):
        """
        Performs a differential expression analysis to identified alterated genes
        and returns a string of the resulting file name
        
        :param fdr_cutoff: the cutoff of FDR to keep the genes
        """
        
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri
        robjects.r('library("edgeR")')
        out_file = "{}.edgeR.txt".format(self.data)
        count = self.data.T
        groups = self.label.unique()
        dgelist = robjects.r.DGEList(counts = count
                      , group = self.label
                      , genes = self.genes
        )

        norm = robjects.r.calcNormFactors(dgelist)
        comdis = robjects.r.estimateCommonDisp(norm)
        y = robjects.r.estimateTagwiseDisp(comdis)
        et = robjects.r.exactTest(y, pair=groups)
        robjects.r.write.table(robjects.r.topTags(et, n = "all"), p.value=fdr_cutoff, out_file, sep="\t")
        return(out_file)

    def testES(self, es=0.8, method='wilcoxon', input_file):
        pass

    # to be implemented
    def mRMR(self, input_file):
        pass

    def fsRF(self, n_trees=10000, input_file):
        """
        Performs feature importance evaluation using random forest
        and returns a list contains the importance value of each gene
        
        :param n_trees: the number of generated trees
        :param input_file: an input data table file with gene expression level 
        """
        
        estimator = RandomForestClassifier(n_estimators=n_trees, random_state=0)

        estimator.fit(self.data, self.label)
        importances = zip(self.genes, estimator.feature_importances_)
        return(sorted(importances, key=lambda x: x[1], reverse=True))

    def fsDT(self, n_trees=10000, input_file):
        """
        Performs feature importance evaluation using randomized decision trees
        and returns a list contains the importance value of each gene
        
        :param n_trees: the number of generated trees
        :param input_file: an input data table file with gene expression level
        """
        
        estimator = ExtraTreesClassifier(n_estimators=25000, random_state=0)
        estimator.fit(self.data, self.label)
        importances = zip(self.genes, estimator.feature_importances_)
        return(sorted(importances, key=lambda x: x[1], reverse=True))

    def fsRFE(self, cv=5, input_file):
        """
        Performs feature importance evaluation using SVM with linear kernel
        and returns a list contains the selected genes by SVC
        
        :param cv: the fold of cross-validation
        :param input_file: an input data table file with gene expression level 
        """
        
        from sklearn.svm import SVC
        from sklearn.feature_selection import RFECV
        estimator = SVC(kernel="linear")
        selector = RFECV(estimator, step=1, cv=cv)
        selector.fit(self.data, self.label)
        markers = selector.support_
        return(self.data.columns[markers])

# Classification model on training data using selected genes
# for prediction of test data 
class Classifier:
    """
    Usage::
    >>> clf = Classifier(training, test)

    """

    def __init__(self, training, test):
        self.training = training
        self.test = test
        self.training_label = self.training.index.get_level_values('class')


    def clfRF(self, n_trees=10000):
        rfc = RandomForestClassifier(n_estimators=n_trees)
        rfc.fit(self.training, self.training_label)
        return(rfc.predict(self.test))

    def clfSVM(self, cv=5, c=[2**i for i in arange(-5, 7, 0.1)], g=[0, 10]):
        from sklearn.svm import SVC
        from sklearn.grid_search import GridSearchCV
        tuned_parameters = [
                   {'kernel': ['linear'], 'C': c}
                ]

        scores = ['precision', 'f1', 'accuracy', 'roc_auc']
        for score in scores:
            clf = GridSearchCV(SVC(C=1), tuned_parameters, cv=4, scoring=score)
            clf.fit(self.training, self.training_label)
            print(clf.best_estimator_)


    def clfLDA(self):
        pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("training",
                        help="the input expression table for \
                        training dataset with count or abundance of genes",
                        type=str)
    parser.add_argument("test",
                        help="the input expression table for \
                        test dataset with count or abundance of genes",
                        type=str)
    parser.add_argument("-s", "--subclass", action='store_true',
                        help="with subclass",
                        default=False
                        )

    args = parser.parse_args()

    training_file = args.training
    test_file = args.test
    training_dataset = dp.Data(training_file, subclass=args.subclass).df
    test_dataset = dp.Data(test_file, dataset='test', subclass=args.subclass).df
    training_label = training_dataset.index.get_level_values('class')
    test_label = test_dataset.index.get_level_values('class')

    frames = (training_dataset, test_dataset)
    whole_dataset = dp.pd.concat(frames)
    whole_stand = dp.standardization(whole_dataset, method='max_min')

    training_stand = whole_stand.xs(
                        'training',
                        level='dataset',
                        drop_level=False
                    )
    test_stand = whole_stand.xs('test', level='dataset', drop_level=False)

    fs1 = FS(training_stand, training_label)
    deg_file = fs1.edgeR(fs1)

    # furthur code for customized gene selection pipeline

if __name__ == '__main__':
    main()
