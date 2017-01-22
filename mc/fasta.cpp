/*
2015-10-30
Updated by Kristin to include more analysis.

2015-11-05
Updated by Kristin to include possibility to have Gouy-Chapman surface and sticky wall potential. Note that the distance between the protein and the surface is calculated differently from before since the surface has been moved in the python script.

2016-01-21
Analysis of adsorbed proteins and protein end-to-end-distances added by Kristin.

2016-06-23
Changed definition of adsorbed proteins to < 20 AA from the surface (for calculation of xy-distance between adsorbed proteins).
Samples the end-to-end distances of all proteins individually.

*/

#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboidslit> Tspace;
//typedef Space<Geometry::Sphere> Tspace;

// Surface potentials
typedef Potential::GouyChapman<> Textpot1;
typedef Potential::StickyWall<> Textpot2;
// Other potentials
typedef Potential::DebyeHuckel Tppel;
typedef Potential::LennardJones Tppsr;
typedef Potential::CombinedPairPotential<Tppel, Tppsr> Tpairpot;

int main()
{

    InputMap mcp("fasta.json");        // Open input parameter file

    Tspace spc(mcp);                   // Simulation space (all particles and group info)

    FormatPQR::save("init.pqr", spc.p, spc.geo.len);

    auto pot = Energy::Nonbonded<Tspace, Tpairpot>(mcp)
        + Energy::Bonded<Tspace>()
        + Energy::EquilibriumEnergy<Tspace>(mcp)
        + Energy::ExternalPotential<Tspace, Textpot1>(mcp)
        + Energy::ExternalPotential<Tspace, Textpot2>(mcp);
    pot.first.second.expot.setSurfPositionZ(&spc.geo.len_half.z());
    pot.second.expot.setSurfPositionZ(&spc.geo.len_half.z());

    /*
     * At this stage all bonds are already setup. We now replace the bond
     * potential with a harmonic potential to which we subtract the short-range
     * part of the nonbonded potential. The net effect is that bonded pairs
     * interact with harmonic+debyehuckel.
     */
    auto bonded = std::get<1>(pot.tuple());
    auto b = bonded->getBondList(); // copy bond info
    bonded->clear();                // clear bonded energy
    for ( auto &i : b )
    {             // and rebuild...
        bonded->add(i.first.first, i.first.second,
                    Potential::Harmonic(mcp["energy"]["nonbonded"]) - 0.9999 * Tppsr(mcp["energy"]["nonbonded"]));
    }

    spc.load("state");                                  // Load start configuration, if any

    Move::Propagator<Tspace> mv(mcp, pot, spc);           // assembly of all MC moves
    Analysis::CombinedAnalysis analysis(mcp, pot, spc);
    Analysis::RadialDistribution<> rdf_ads_prot(0.4);   // adsorbed protein-protein rdf (xy)
    Histogram<double> rdf(0.4);                         // protein-surface rdf
    std::map<int, Analysis::LineDistribution<> > surfresdist; // residue-surface rdf
    Table2D<double, Average<double> > netqtableProtein(0.4); // protein charge
    Table2D<double, Average<double> > netqtableSurface(0.4); // surface charge
    Table2D<double, Average<double> > rg2table(0.4);         // Rg
    Table2D<double, Average<double> > rgx2table(0.4);
    Table2D<double, Average<double> > rgy2table(0.4);
    Table2D<double, Average<double> > rgz2table(0.4);

    cout << spc.info() + pot.info() + textio::header("MC Simulation Begins!");

    auto pol = spc.findMolecules("protein");
    auto sur = spc.findMolecules("surface");
    Table2D<double, Average<double> > ree2table[pol.size()];
    int scnt = 0;

    auto shape = analysis.get<Analysis::PolymerShape<Tspace>>();
    if (shape==nullptr)
        throw std::runtime_error("Add PolymerShape to analysis via JSON input");

    MCLoop loop(mcp);    // handle mc loops
    while ( loop[0] )
    {
        while ( loop[1] )
        {
            mv.move();
            analysis.sample();

            if ( slump() > 0.95 )
            {

                int prot_index = 0; // The index of the protein that we are looking at in the loop.
                for ( Group* i : pol )
                {

                    double distance = spc.geo.len.z() / 2 - i->cm.z(); // distance between protein and surface
                    rdf(distance)++;                             // sample protein-surface distribution
                    netqtableProtein(distance) += netCharge(spc.p, *i);
                    Point rg2 = shape->vectorgyrationRadiusSquared(*i, spc);
                    rg2table(distance) += rg2.x() + rg2.y() + rg2.z();
                    rgx2table(distance) += rg2.x();
                    rgy2table(distance) += rg2.y();
                    rgz2table(distance) += rg2.z();

                    // Track end-to-end distances of the individual proteins.
                    double ree2 = spc.geo.sqdist(spc.p[i->front()], spc.p[i->back()]);
                    ree2table[prot_index](distance) += ree2;

                    // Sample protein-protein distances (xy-direction) on surface.
                    if ( distance < 20 )
                    {
                        Point position_i(i->cm.x(), i->cm.y(), 0);
                        for ( Group* i2 : pol )
                        {
                            if ((spc.geo.len.z() / 2 - i2->cm.z()) < 20 && i != i2 )
                            {
                                Point position_i2(i2->cm.x(), i2->cm.y(), 0);
                                rdf_ads_prot(spc.geo.dist(position_i, position_i2))++; // double counting.
                            }
                        }
                    }

                    // Sample residue-surface distribution. Should work for several proteins.
                    for ( int j : *i )
                    {
                        double resdist = spc.geo.len.z() / 2 - spc.p[j].z();
                        surfresdist[j % i->size()](resdist)++;
                    }
                    prot_index++;
                }
                for ( auto i : sur )
                {
                    for ( Group* j : pol )
                    {
                        double distance = spc.geo.len.z() / 2 - j->cm.z();
                        // Sample net charge of surface as function of distance to protein. Probably only makes sense with one protein.
                        netqtableSurface(distance) += ::netCharge(spc.p, *i);
                    }
                }
                scnt++;
            }

        } // end of micro loop

        cout << loop.timing();                                 // print timing and ETA information

    } // end of macro loop

    netqtableProtein.save("netq_protein_dist.dat");
    netqtableSurface.save("netq_surface_dist.dat");
    rg2table.save("rg2_dist.dat");
    rgx2table.save("rgx2_dist.dat");
    rgy2table.save("rgy2_dist.dat");
    rgz2table.save("rgz2_dist.dat");
    rdf.save("rdf.dat");
    rdf_ads_prot.save("rdf_ads_prot.dat");
    for ( unsigned int i = 0; i < pol.size(); i++ )
    {
        std::string filename = "ree2_prot" + std::to_string(i + 1) + "_dist.dat";
        ree2table[i].save(filename);
    }

    std::ofstream file("surf_res_dist.dat");
    for ( double d = 0; d <= spc.geo.len.z(); d += 0.2 )
    {
        for ( int i = 0; i < pol.front()->size(); i++ )
            file << d << " " << i + 1 << " " << -log(double(surfresdist[i](d)) / (double(scnt) * pol.size()))
                 << endl; // Note: This will write 'inf' if probability is 0
        file << endl;
    }
    file.close();

    cout << loop.info() + mv.info() + analysis.info();
}
