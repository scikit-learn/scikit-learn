__all__ = ['tree','algorithm','ga_util','ga_list','gene','genome','language',
           'population','scaling','selection','tree_opt']

for module in __all__:
    exec('import ' + module)           
# need to look at ga_gnm,tree_opt to see if they really fit.