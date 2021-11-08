"""
===========================
Movie Recommendation Engine
===========================
This example shows how scikit-learn can be used to build a
content-based recommendation engine.
"""

# Author: Sven Eschlbeck <sven dot eschlbeck at t-online dot de>
# License: BSD 3 clause

# Standard scientific Python imports
import pandas as pd
import numpy as np

# Import datasets, classifiers and performance metrics
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity

#################################################################
# Recommendation Engines
#-----------------------
# A Recommendation Engine seeks to predict the rating a user would
# give to an item or to predict the item in order to meet the
# user's taste. The more precise Recommendation Engines are, the
# better they can guess a user's or customer's taste.
# Recommendation Engines are used across various scientific and
# economic fields, e.g. in playlist generators, product recommen-
# ders in online stores or content recommenders on social media
# platforms.
#
# Although all Recommendation Engines aim to achieve the same goal,
# they use different techniques. The three main techniques are:
# - Popularity-based
# - Content-based
# - Collaborative Filtering
#
# Popularity-based
#-----------------
# Popularity-based Recommendation Engines are the simplest form of
# Recommendation Engine. They count the number of views or clicks
# a product or service got and rank their popularity alongside
# other products and services. The user then gets recommended the
# currently most popular item(s). An example for Popularity-based
# Recommendation Engines is YouTube's list of trending videos.
#
# Content-based
#--------------
# Content-based Recommendation Engines focus on the preferences
# of the user currently using the service. They take into con-
# sideration what content the user likes, then make out the key
# components or features of the product or service and finally
# search for other products and services with similar characte-
# ristics to be recommended to that user.
#
# Collaborative Filtering
#------------------------
# Collabarative Filtering tries to identify similar types of
# users based on their individual behavior and preferences.
# If two users were clustered as members of a group of similar
# users, the Recommendation Engine can suggest content that one
# group member liked, to the other group members. Since they
# have shown a related consumption behavior befor, chances are
# hight that they will like the recommendation, too.

#################################################################
# Content-based Movie Recommendation Engine
#------------------------------------------
# In this example, our Recommendation Engine will be Content-based.
# Consequently, we need to find similar movies to a given movie
# in order to make a movie proposal the user might like.

# Read .csv file and write data into dataframe
df = pd.read_csv("movies.csv")

# Choose relevant columns of the dataframe as feature set
features = ["keywords", "cast", "genres", "director"]

# Concatenating all extracted features into single string
def concat_features(row):
	return (row["keywords"] + " " + row["cast"] + " " +
		row["genres"] + " " + row["director"])

# Replacing NaN values with blank strings
for feature in features:
	df[feature] = df[feature].fillna("")

# Iterating over each dataframe row, applying concat_features
# method to each row and writing resulting string in the
# newly created concatenated_features column
df["concatenated_features"] = df.apply(concat_features, axis=1)

# Calling CountVectorizer to receive the count matrix
# listing the number of occurences of each word/element
cv = CountVectorizer()
count_matrix = cv.fit_transform(df["concatenated_features"])

# Calculate the cosine similarity matrix from the count matrix
# All word matrices base on the same principle. Word similarity
# is best representable by imagining a multi-dimensional
# vector space. Each phrase or word is represented by a unique
# vector inside the vector space. Similarity between two words
# is then defined by the distance to each other. Close vectors
# stand for related phrases/words. Insted of simply defining
# the distance between two vectors as their angular distance,
# the cosine(theta) with theta being the angle between the
# vectors is used.
cosine_sim = cosine_similarity(count_matrix)

# Retrieving a movie title from an index
def get_title_from_index(index):
	return df[df.index == index]["title"].values[0]

# Retrieving a movie index from a title
def get_index_from_title(title):
	return df[df.title == title]["index"].values[0]

# Ask user for the title of movie he/she likes
movie_user_likes = "Avatar"

# Get index of that movie with pre-defined function
movie_index = get_index_from_title(movie_user_likes)

# Jump to the row belonging to the user movie in the similarity
# matrix. This matrix lists the similarity of that movie to all
# other movies. Enumerate through similarity scores and create
# tuples of form (movie index, similarity score).
similar_movies = list(enumerate(cosine_sim[movie_index]))

# Sort list similar_movies by similarity scores in descending
# order. Drop first element after sorting, since the most
# similar movie to a given movie is itself.
sorted_similar_movies = sorted(similar_movies, key=lambda x:x[1], reverse=True)[1:]

# Print 5 most recommended movies (first 5 list entries)
i=0

for element in sorted_similar_movies:
	print(get_title_from_index(element[0]))
	i+=1
	if i>5:
		break