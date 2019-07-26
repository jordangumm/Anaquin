#ifndef RESOURCES_HPP
#define RESOURCES_HPP

#include "data/data.hpp"

namespace Anaquin
{
    Scripts PlotGene();
    Scripts PlotInsert();
    Scripts PlotLinear();
    Scripts PlotAllele();
    Scripts PlotLDensity();
    Scripts PlotKSomatic();
    Scripts PlotKSynthetic();
    Scripts PlotStrelkaROC();
    Scripts PlotQualFilter();

    Scripts rna();
    Scripts meta();
    Scripts split();
    Scripts somatic();
    Scripts germline();
    Scripts partition();
    Scripts calibrate();
    
    Scripts CSS();
    Scripts PReport();
    Scripts GermHTML();
    Scripts SomaHTML();
    Scripts GSplitHTML();
    Scripts RSplitHTML();
    Scripts MSplitHTML();
    Scripts CalibrateHTML();
}

#endif
