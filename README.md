This tutorial was developed for the SNP Rotation Seminar given by Adam Leaché, University of Washington, Winter 2022. 

The purpose is to demonstrate DAPC and ADMIXTURE visualizations in R. 

# Background

The neotropical tree species *Cedrela odorata* is a target of illegal logging and has been documented in regional timber trade for over 250 years (Pennington and Muellner, 2010). Of these species, *Cedrela odorata* (Spanish cedar; cedro), is the second most traded neotropical tree species, following bigleaf mahogany (*Swietenia macrophylla*). In 2001, *C. odorata* was listed under the protections of CITES Appendix III requiring validated documentation of species identity and source for both export and import documentation, protecting populations in Bolivia, Brazil, Colombia, Guatemala, and Peru (Ferriss, 2014). Consequently, two other species (*C. fissilis* Vell. and *C. lilloi syn. angustifolia*) have also been afforded protection under CITES Appendix III, since their survival is threatened due to strong similarities with *C. odorata* and its trade-demand (Ferriss, 2014; UNEP-WCMC, 2015). Despite the anticipated increase in protection for *C. odorata* and its “look-alikes”, law enforcement of CITES regulations and the continuance of CITES protections requires credible estimates of species presence in trade, accurate geographic distributions, and up-to-date extinction risk (Text of the Convention; Contracting States, 1973). However, differentiating *Cedrela* species is challenging, and differentiating species based on wood tissue alone is nearly impossible, even for experts. The data presented here are part of a larger project to provide genomic resources for *Cedrela*, including a SNP genotyping assay for species and origin identification.

Although the geographic distribution of *Cedrela odorata* is widespread, occurring from Mexico to northern Argentina (24º N to 27º S) and including the islands of the Caribbean (Pennington and Muellner, 2010), it primarily exists in low-density stands with maybe one individual per hectare in Peru, Costa Rica, Colombia, and Guyana (Tajikistan 2019). Strong light dependency and pests in its native range lead to poor natural regeneration of stands and fragmented populations (Cintron, 1990; Pennington and Muellner, 2010). For example, reduced population connectivity may limit the spread of regional genotypes that are tolerant to the mahogany shoot borer, *Hypsipyla grandella* Zeller (Newton et al., 1999). Deforestation has led to a nearly 30% decline in the global distribution of *C. odorata* in the last 100 years, and that estimate is expected to increase by over 10% in the next 100 years (Tajikistan 2019).The primary threat to *C. odorata* and other *Cedrela* is decreased population connectivity from exploitative logging (Muellner et al., 2010, 2011; Pennington and Muellner, 2010; Cavers et al., 2013), and the demand for *Cedrela* wood is insatible.

Oddly, *Cedrela* grow very well outside of their native range, even being considered invasive in some areas. This could be a example of the, "enemy release hypothesis," in action, where species escape the primary threat to their population growth outside of their native range (e.g., the mahogany shoot borer). One place where *Cedrela* are considered invasive is the Galapagos archipeligo (Ecuador). Part of this project is to determine the source of this population, and this is the data presented here and the reseach question for this session. 

Data presented here are next generation sequence data from dissertation research (KNF advised by: Andy Jones [Oregon State University], Rich Cronn [USDA Forest Service PNW Research Station]) and an active collaboration with
 Gonzalo Rivas-Torres and Maria de Lordes Torres of Universidad San Francisco de Quito. Dissertation specimens were collected from the Missouri Botanical Garden Herbarium, and Galapagos/Ecuador specimens of interest were provided by Rivas-Torres and de Lourdes-Torres. Specimen sequence data was either generated by target capture, amplicon sequencing, or both.

The full specimen alignment contained 282 specimens and nearly 1 million SNPs. Variant call format file (vcf) was generated with BCFtools mpileup Filtering was performed with VCFtools: (biallelic, minimum depth = 100, 0% missing data, and 1 SNP per 10kb). For the exercise, I focused on Ecuador. 

The vcf provided here contains 127 individuals and 725 SNPs. 

# Files included

File | Description
-----|------------
20220213_snp_seminar_fil_thin.recode.vcf | VCF file (127 individuals and 725 SNPs)
20220213_snp_seminar_samples.csv | meta data for samples
country_labels_2020.csv | country labels for maps
map_files | directory containing shapefiles from http://thematicmapping.org/
snp_seminar.R | the R script of this tutorial
snp_seminar.html |the R markdown of this tutorial rendered as an HTML - download to view
snp_seminar_admixture |directory containing ADMIXTURE Outputs
snp_seminar_packages.R | list of packages and code to install and load them

# Literature Cited

Alexander, D. H., J. Novembre, and K. Lange. 2009. Fast model-based estimation of ancestry in unrelated individuals. Genome Research 19: 1655–1664.

Cintron, B. B. 1990. *Cedrela odorata* L. Cedro hembra, Spanish cedar. Silvics of North America 2: 250.

Ferriss, S. 2014. An analysis of trade in five CITES-listed taxa. The Royal Institute of International Affairs, Chatham House, London, England, UK.

Newton, A. C., T. R. Allnutt, A. C. M. Gillies, A. J. Lowe, R. A. Ennos, A. C. Newton, T. R. Allnutt, et al. 1999. Molecular phylogeography, intraspecific variation and the conservation of tree species. Trends in Ecology & Evolution 14: 140–145.

Pennington, T. D., and A. N. Muellner. 2010. A monograph of *Cedrela* (Meliaceae). dh books, Milborne Port, England, UK.

Tajikistan, M. 2019. Proposal for amendment of Appendix I or II for CITES CoP18 Prop. 57. Consideration of proposals for amendment of appendices I and II, 1–26. CITES, Colombo, Sri Lanka.

UNEP-WCMC. 2015. Overview of CITES Appendix III listings. UNEP-WCMC, Cambridge, U. K.
