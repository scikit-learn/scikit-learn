<!--
Thanks for contributing a pull request! Please ensure you have taken a look at
the contribution guidelines: https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md#pull-request-checklist
-->

#### Reference Issues/PRs

Issue Number: 13990
Issue: Add 3-fold split method train/val/test

Issue Link: https://github.com/scikit-learn/scikit-learn/issues/13990
We had created issue number 16342 which was closed and marked as duplicate
Link For Closed Issue:https://github.com/scikit-learn/scikit-learn/issues/16342 


#### What does this implement/fix? Explain your changes.

We have introduced an additional functionality in the existing function of train_test_split. Now the Target(Y) and Response(X) can be split three ways i.e. into train set / validation set / test set. 
This is an additional optional functionality and does not hamper the existing functionality.  
	- only if function is called with the optional parameter "validation_size", then the function returns 3 subsets. Otherwise the existing functionality to return two subsets for each input array is retained. 

The test_split.py is also changed in order to test the new functionality 

#### Any other comments?

Fixes were made in the file fixes.py to resolve import issues due to different scipy versions. 
Source Code link: https://github.com/rbewoor/scikit-learn
Different use cases were tested using the script Test_Cases_SciKit-Learn.ipynb in the folder Dedicated_Test_Cases
Link https://github.com/rbewoor/scikit-learn/blob/master/Dedicated_Test_Cases/Test_Cases_SciKit-Learn.ipynb
