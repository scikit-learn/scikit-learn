from scipy import ga

def fitness_ex1(self):
    score = 0
    desired = ['a','b','c','d','b',-1,.25,.6, .8, 0]
    # compare letters
    actual = self.get_values()
    pairs = zip(actual[:5],desired[:5])
    for g,d in pairs:
        if g == d:
            score += 1
    pairs = zip(actual[5:],desired[5:])
    for g,d in pairs:
        diff = g - d
        score -= abs(diff)
    return score        
    
def ex1():
    """ This example illustrates the main steps in setting
        up a genetic optimization:
        
          1. Specify the genes types used to encode your problem
          2. Group these genes into a genome.
             a.  Specify the fitness function that evaluates the genomes.
          3. Create a population of the genomes
          4. Specify the algorithm used to evolve the population
          5.  
    """
    
    # 1. First scpecify your genes.  To gene types are 
    #    currently supported, list_gene and float_gene
    
    #    A list gene chooses its value from a list of values.
    #    The list can contain any type of object.
    g1 = ga.gene.list_gene(['a','b','c','d'])    
    #    Float genes take on a continuous value between two
    #    bounds.
    g2 = ga.gene.float_gene((-1.,1.))    
    #    We'll replicate these genes several times to make a longer
    #    genome.
    all_genes = g1.replicate(5) + g2.replicate(5)    
    # 2.  Create a specialized "list_genome" (as opposed to tree_genome) 
    #     class with the desired fitness function.
    #     It's structure is defined by our gene list. 
    class this_genome(ga.genome.list_genome):
        pass
    this_genome.performance = fitness_ex1            
    gnm = this_genome(all_genes)

    # 3.  Create a population of the genomes.
    #
    pop = ga.population.population(gnm)
    # 4.  Now use the basic genetic algorithm to evolve the population
    #
    galg = ga.algorithm.galg(pop)
    # change a few settings
    settings = {'pop_size':250,'p_replace':.8,'p_cross': .8, 'p_mutate':'gene',
                'p_deviation': 0.,'gens':35,'rand_seed':0,'rand_alg':'CMRG'}
    galg.settings.update(settings)
    galg.evolve()
    print galg.pop.best()
    
if __name__ == '__main__':
    ex1()    