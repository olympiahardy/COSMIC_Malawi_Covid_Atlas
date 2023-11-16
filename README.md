# Covid Study using Minimally invasive autopsy for Interrogation of Cellular response (COSMIC) in Malawi to generate a spatially resolved single cell Covid atlas
Codebase for a spatially resolved single-cell atlas of the lung in fatal Covid19 in an African population that reveals a distinct cellular signature and an interferon gamma dominated response.

James Nyirenda*, Olympia Hardy*, João Da Silva Filho*, Vanessa Herder, Charalampos Attipa, Charles Ndovi, Memory Siwombo, Takondwa Namalima, Leticia Suwedi, Watipenge Nyasulu, Thokozile Ngulube, Deborah Nyirenda, Leonard Mvaya, Joseph Phiri, Dennis Chasweka, Chisomo Eneya, Chikondi Makwinja, Chisomo Phiri, Frank Ziwoya, Abel Tembo, Kingsley Makwangwala, Stanley Khoswe, Peter Banda, Ben Morton, Orla Hilton, Sarah Lawrence, Monique Freire dos Reis, Gisely Cardoso Melo, Marcus Vinicius Guimaraes de Lacerda, Fabio Trindade Maranhão Costa, Wuelton Marcelo Monteiro, Luiz Carlos de Lima Ferreira, Carla Johnson, Dagmara McGuinness, Kondwani Jambo, Michael Haley, Benjamin Kumwenda, Massimo Palmarini, Kayla G. Barnes+, Donna M. Denno+, Wieger Voskuijl+ , Steve Kamiza+, Kevin Couper+, Matthias Marti+, Thomas Otto+, Christopher A. Moxon+,**

# Check out our pre-print on Biorxiv:
TBA

# Background
Postmortem single-cell studies have transformed understanding of lower respiratory tract diseases (LRTD) including Covid19 but there is almost no data from African settings where HIV, malaria and other environmental exposures may affect disease pathobiology and treatment targets. We used histology and high-dimensional imaging to characterise fatal lung disease in Malawian adults with (n=9) and without (n=7) Covid19, and generated single-cell transcriptomics data from lung, blood and nasal cells. Data integration with other cohorts showed a conserved Covid19 histopathological signature, driven by contrasting immune and inflammatory mechanisms: in the Malawi cohort, by response to interferon-gamma (IFN-y) in lung-resident alveolar macrophages, in USA, European and Asian cohorts by type I/III interferon responses, particularly in blood-derived monocytes. HIV status had minimal impact on histology or immunopathology. Our study provides data resources and highlights the importance of studying the cellular mechanisms of disease in underrepresented populations, indicating shared and distinct targets for treatment.

# Atlas Resources
Check out our interactive single cell atlases hosted on the cellxgeneVIP platform hosted by the University of Glasgow, by clicking on the URLs below:
* Lung Atlas - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_Lung_Atlas.h5ad/
* Lung Immune Atlas - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_Lung_Immune_Atlas.h5ad/
* Lung Stromal Atlas - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_Lung_Stromal_Atlas.h5ad/
* Nasal Atlas - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_Nasal_Atlas.h5ad/
* Blood Atlas - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_Blood_Atlas.h5ad/
* Histopathology slides on virtual microscope: https://covid-atlas.cvr.gla.ac.uk
* IMC - https://cellatlas-cxg.mvls.gla.ac.uk/COSMIC/view/COSMIC_IMC_Lung.h5ad/

# Code Layout
  * Analysis
      * scRNA - Scripts containing all pre-processing steps to generate single cell RNA sequencing objects for each tissue
      * IMC - Scripts containing all pre-processing steps to generate the imaging mass cytometry object for the lung
  * Figures
      * Scripts containing all steps to regenerate figures from the study
   
# Raw and processed data access:
TBA (EBI ArrayExpress)
TBA (Zenodo)
