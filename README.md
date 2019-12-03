# OpenBreast
This is a Matlab implementation of different tools for processing digital mammography images developed by Universidad Industrial de Santander. OpenBreast was publicly released in [1] and has been clinically evaluated for the task of breast cancer risk assessment in [2]. The following tasks have been implemented:

* Feature extraction for parenchymal analysis [1]
* Image standardization for (RAW and PROCESSED) digital mammography images
* Breast segmentation and chest wall detection [3]
* Detection of regions on interest within the breast [4,5]
* Breast density segmentation [6]

To get started first run setup.m to configure Openbreast. Then run the following demos:

* demo01		Breast segmentation
* demo02 		ST mapping
* demo03 		ROI detection
* demo04 		Feature extraction on FFDM images
* demo05 		Breast density segmentation

For further details, please refer to: https://sites.google.com/view/cvia/openbreast. OpenBreast is also available at FileExchange [![View OpenBreast on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/71615-openbreast)

[1] S. Pertuz et al., Open Framework for Mammography-based Breast Cancer Risk Assessment, IEEE-EMBS International Conference on Biomedical and Health Informatics, 2019.

[2] S. Pertuz et al., Clinical evaluation of a fully-automated parenchymal analysis software for breast cancer risk assessment: A pilot study in a Finnish sample,
European Journal of Radiology: 121, 2019.

[3] B. Keller et al., Estimation of breast percent density in raw and processed full field digital mammography images via adaptive fuzzy c-means clustering and support vector machine segmentation, Med. Phys, 2012.

[4] S. Pertuz, C. Julia, D. Puig, A novel mammography image representation framework with application to image registration, Proc. International Conference on Pattern Recognition, 2014.

[5] G. Torres, S. Pertuz, Automatic Detection of the Retroareolar Region in Mammograms, Proc. Latin American Congress on Biomedical Engineering, 2016

[6] G. F. Torres et al., "Morphological Area Gradient: System-independent Dense Tissue Segmentation in Mammography Images," Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC), 2019.


