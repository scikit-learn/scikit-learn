#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import pandas as pd

"""
    - COPYRIGHT: Thomas E. Cason, Undergraduate Research Assistant
        Division of Biostatistics and Epidemiology
        Department of Health Evaluation Sciences
        University of Virginia School of Medicine
        Box 600, Charlottesville, VA 22908 USA
        Electronic Mail:  tcason@virginia.edu
    - SOURCE: Hind, Philip.  "Encyclopedia Titanica."  Online.  Internet.
    n.p.  02 Aug 1999.  Avaliable http://atschool.eduweb.co.uk/phind
    http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3.csv
    - DESCSHORT: The titanic3 data frame describes the survival
status of individual passengers on the Titanic.  The titanic3 data
frame does not contain information for the crew, but it does contain
actual and estimated ages for almost 80%% of the passengers.
    - DESCLONG: http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3info.txt
    - NOTE:
        pclass          Passenger Class
                        (1 = 1st; 2 = 2nd; 3 = 3rd)
        survival        Survival
                        (0 = No; 1 = Yes)
        name            Name
        sex             Sex
        age             Age
        sibsp           Number of Siblings/Spouses Aboard
        parch           Number of Parents/Children Aboard
        ticket          Ticket Number
        fare            Passenger Fare
        cabin           Cabin
        embarked        Port of Embarkation
                        (C = Cherbourg; Q = Queenstown; S = Southampton)
        boat            Lifeboat
        body            Body Identification Number
        home.dest       Home/Destination

    - SPECIAL NOTES:
        Pclass is a proxy for socio-economic status (SES)
         1st ~ Upper; 2nd ~ Middle; 3rd ~ Lower

        Age is in Years; Fractional if Age less than One (1)
         If the Age is Estimated, it is in the form xx.5

        Fare is in Pre-1970 British Pounds (£)
         Conversion Factors: £1 = 12s = 240d and 1s = 20d

        With respect to the family relation variables (i.e. sibsp and parch)
        some relations were ignored.  The following are the definitions used
        for sibsp and parch.

        Sibling:  Brother, Sister, Stepbrother, or Stepsister of Passenger Aboard Titanic
        Spouse:   Husband or Wife of Passenger Aboard Titanic (Mistresses and Fiancées Ignored)
        Parent:   Mother or Father of Passenger Aboard Titanic
        Child:    Son, Daughter, Stepson, or Stepdaughter of Passenger Aboard Titanic

        Other family relatives excluded from this study include cousins,
        nephews/nieces, aunts/uncles, and in-laws.  Some children travelled
        only with a nanny, therefore parch=0 for them.  As well, some
        travelled with very close friends or neighbors in a village, however,
        the definitions do not support such relations.

"""


def load():
    """Returns a pandas DataFrame containing the data from titanic3"""
    return pd.read_csv("%s/titanic3.csv" %
                       os.path.dirname(os.path.abspath(__file__)),
                       na_values=None)
