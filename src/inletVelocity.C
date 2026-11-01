// OpenFOAM codedFixedValue body, used on both freestream boundaries.
const auto c = diamond::conditions();
operator==(Foam::vector(c.ux, c.uy, 0));
