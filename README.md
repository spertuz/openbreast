# OpenBreast
This is a Matlab implementation of different tools for processing digital mammography images developed by Universidad Industrial de Santander. The following tasks have been implemented [1]:

* Image standardization for (RAW and PROCESSED) digital mammography images
* Breast segmentation and chest wall detection [2]
* Detection of regions on interest within the breast [3,4]
* Feature extraction for parenchymal analysis [1]

To get started first run setup.m to configure Openbreast. Then run the following demos:

* demo01		Breast segmentation
* demo02 		ST mapping
* demo03 		ROI detection
* demo04 		Feature extraction on FFDM images

For further details, please refer to: https://sites.google.com/view/cvia/openbreast

[1] S. Pertuz et al., Open Framework for Mammography-based Breast Cancer Risk Assessment, IEEE-EMBS International Conference on Biomedical and Health Informatics, 2019.
[2] B. Keller et al., Estimation of breast percent density in raw and processed full field digital mammography images via adaptive fuzzy c-means clustering and support vector machine segmentation, Med. Phys, 2012.
[3] S. Pertuz, C. Julia, D. Puig, A novel mammography image representation framework with application to image registration, Proc. International Conference on Pattern Recognition, 2014.
[4] G. Torres, S. Pertuz, Automatic Detection of the Retroareolar Region in Mammograms, Proc. Latin American Congress on Biomedical Engineering, 2016
