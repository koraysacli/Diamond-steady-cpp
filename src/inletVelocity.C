// OpenFOAM codedFixedValue govdesi. Her iki serbest akis sinirinda kullanilir.
const auto c = diamond::conditions();
operator==(Foam::vector(c.ux, c.uy, 0));
