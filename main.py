from predict import *
from similarity import *
from validate import *
import os
import logging
logging.disable()


def main():
    # create results folder
    os.system("mkdir results")
    os.system("mkdir results/predicted_scores")
    os.system("mkdir results/similarity_scores")
    os.system("mkdir results/gaps")

    # predict scores for reactions in reaction pool
    repeat = 5
    for i in range(repeat):
        get_prediction_score(name='zimmermann')

    # predict mean similarity between candidate reactions and existing reactions
    get_similarity_score(name='zimmermann', top_N=2000)

    # predict metabolic phenotypes
    # If you only want prediction and similarity scores, comment out the following line
    validate()


if __name__ == "__main__":
    main()
