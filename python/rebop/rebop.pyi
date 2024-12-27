import xarray

class Gillespie:
    """Reaction system composed of species and reactions."""

    def add_reaction(
        self,
        /,
        rate: float,
        reactants: list[str],
        products: list[str],
        reverse_rate: float | None = None,
    ) -> None:
        """Add a Law of Mass Action reaction to the system.

        The forward reaction rate is `rate`, while `reactants` and `products` are
        lists of respectively reactant names and product names.
        Add the reverse reaction with the rate `reverse_rate` if it is not `None`.
        """

    def nb_reactions(self, /) -> int:
        """Number of reactions currently in the system."""

    def nb_species(self, /) -> int:
        """Number of species currently in the system."""

    def run(
        self,
        init: dict[str, int],
        tmax: float,
        nb_steps: int,
        seed: int | None = None,
    ) -> xarray.Dataset:
        """Run the system until `tmax` with `nb_steps` steps.

        The initial configuration is specified in the dictionary `init`.
        Returns an xarray Dataset.
        """
