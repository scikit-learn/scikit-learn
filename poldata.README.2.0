
=======

Introduction

This README v2.0 (June, 2004) for the v2.0 polarity dataset comes from
the URL http://www.cs.cornell.edu/people/pabo/movie-review-data .

=======

What's New -- June, 2004

This dataset represents an enhancement of the review corpus v1.0
described in README v1.1: it contains more reviews, and labels were
created with an improved rating-extraction system.

=======

Citation Info 

This data was first used in Bo Pang and Lillian Lee,
``A Sentimental Education: Sentiment Analysis Using Subjectivity Summarization 
Based on Minimum Cuts'',  Proceedings of the ACL, 2004.

@InProceedings{Pang+Lee:04a,
  author =       {Bo Pang and Lillian Lee},
  title =        {A Sentimental Education: Sentiment Analysis Using Subjectivity Summarization Based on Minimum Cuts},
  booktitle =    "Proceedings of the ACL",
  year =         2004
}

=======

Data Format Summary 

- review_polarity.tar.gz: contains this readme and  data used in
  the experiments described in Pang/Lee ACL 2004.

  Specifically:

  Within the folder "txt_sentoken" are the 2000 processed down-cased
  text files used in Pang/Lee ACL 2004; the names of the two
  subdirectories in that folder, "pos" and "neg", indicate the true
  classification (sentiment) of the component files according to our
  automatic rating classifier (see section "Rating Decision" below).

  File names consist of a cross-validation tag plus the name of the
  original html file.  The ten folds used in the Pang/Lee ACL 2004 paper's
  experiments were:

     fold 1: files tagged cv000 through cv099, in numerical order
     fold 2: files tagged cv100 through cv199, in numerical order     
     ...
     fold 10: files tagged cv900 through cv999, in numerical order

  Hence, the file neg/cv114_19501.txt, for example, was labeled as
  negative, served as a member of fold 2, and was extracted from the
  file 19501.html in polarity_html.zip (see below).

  Each line in each text file corresponds to a single sentence, as
  determined by Adwait Ratnaparkhi's sentence boundary detector
  MXTERMINATOR.
 
  Preliminary steps were taken to remove rating information from the
  text files, but only the rating information upon which the rating
  decision was based is guaranteed to have been removed. Thus, if the
  original review contains several instances of rating information,
  potentially given in different forms, those not recognized as valid
  ratings remain part of the review text.
	
- polarity_html.zip: The original source files from which the
  processed, labeled, and (randomly) selected data in
  review_polarity.tar.gz was derived.

  Specifically:  

  This data consists of unprocessed, unlabeled html files from the
  IMDb archive of the rec.arts.movies.reviews newsgroup,
  http://reviews.imdb.com/Reviews. The files in review_polarity.tar.gz
  represent a processed subset of these files. 

=======

Rating Decision (Appendix A)

This section describes how we determined whether a review was positive
or negative.

The original html files do not have consistent formats -- a review may
not have the author's rating with it, and when it does, the rating can
appear at different places in the file in different forms.  We only
recognize some of the more explicit ratings, which are extracted via a
set of ad-hoc rules.  In essence, a file's classification is determined
based on the first rating we were able to identify.


- In order to obtain more accurate rating decisions, the maximum
	rating must be specified explicitly, both for numerical ratings
	and star ratings.  ("8/10", "four out of five", and "OUT OF
	****: ***" are examples of rating indications we recognize.)

- With a five-star system (or compatible number systems):
	three-and-a-half stars and up are considered positive, 
	two stars and below are considered negative.
- With a four-star system (or compatible number system):
	three stars and up are considered positive, 
	one-and-a-half stars and below are considered negative.  
- With a letter grade system:
	B or above is considered positive,
	C- or below is considered negative.

We attempted to recognize half stars, but they are specified in an
especially free way, which makes them difficult to recognize.  Hence,
we may lose a half star very occasionally; but this only results in 2.5
stars in five star system being categorized as negative, which is 
still reasonable.


