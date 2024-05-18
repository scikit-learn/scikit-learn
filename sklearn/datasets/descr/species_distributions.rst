.. _species_distribution_dataset:

Species distribution dataset
----------------------------

This dataset represents the geographic distribution of two species in Central and
South America. The two species are:

- `"Bradypus variegatus" <http://www.iucnredlist.org/details/3038/0>`_ ,
  the Brown-throated Sloth.

 - `"Microryzomys minutus" <http://www.iucnredlist.org/details/13408/0>`_ ,
   also known as the Forest Small Rice Rat, a rodent that lives in Peru,
   Colombia, Ecuador, Peru, and Venezuela.

The dataset is not a typical dataset since a :class:`~sklearn.datasets.base.Bunch`
containing the attributes `data` and `target` is not returned. Instead, we have
information allowing to create a "density" map of the different species.

The grid for the map can be built using the attributes `x_left_lower_corner`,
`y_left_lower_corner`, `Nx`, `Ny` and `grid_size`, which respectively correspond
to the x and y coordinates of the lower left corner of the grid, the number of
points along the x- and y-axis and the size of the step on the grid.

The density at each location of the grid is contained in the `coverage` attribute.

Finally, the `train` and `test` attributes contain information regarding the location
of a species at a specific location.

The dataset is provided by Phillips et. al. (2006).

.. topic:: References

 * `"Maximum entropy modeling of species geographic distributions"
   <http://rob.schapire.net/papers/ecolmod.pdf>`_ S. J. Phillips,
   R. P. Anderson, R. E. Schapire - Ecological Modelling, 190:231-259, 2006.
