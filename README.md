# Ruling Laboratory — private prototype

This is a static, unpublished prototype in its own private repository.  It is
not part of Mingyuan Hu's public homepage; the homepage may link to it only
after the tool has been validated and intentionally published.

The first version accepts a lattice polygon as a coordinate list and recognizes
two fully verified examples:

- Ding–Wei Fig. 2(c): `[(0,0),(4,0),(1,2)]`
- Ding–Wei Fig. 3(p): `[(0,0),(3,0),(3,1),(0,2)]`

For a recognized polygon it shows:

- direct all-disk, annular, and total rational ruling counts;
- the standard ruling;
- rational ruling representatives modulo the selected diagram symmetry;
- multiplicities and the `D/R/A/B/V` profiles.

The data and SVG diagrams are generated from direct ruling certificates in the
`Nodal_Curves` project. Tropical counting is not imported by the exporter.

This repository has not been deployed.
