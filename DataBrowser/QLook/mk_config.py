import json

tab1_dspec_wdth = 1450
tab1_dspec_hght = 500
tab1_dspec_xPro_wdth = tab1_dspec_wdth
tab1_dspec_xPro_hght = 150
tab1_dspec_yPro_wdth = 150
tab1_dspec_yPro_hght = tab1_dspec_hght
tab1_StrID_DataTb_BUT_wdth = 200
tab1_dspec_Ctrl_widget_wdth = 200
tab1_StrID_DataTb_wdth = tab1_dspec_wdth - tab1_StrID_DataTb_BUT_wdth - tab1_dspec_Ctrl_widget_wdth  # - 450
# tab1_StrID_DataTb2_wdth = 450
tab1_StrID_DataTb_hght = 280
# tab1_StrID_DataTb2_hght = tab1_StrID_DataTb_hght

tab2_dspec_wdth = 1290
tab2_dspec_hght = 220
tab2_dspec_xPro_wdth = 900
tab2_dspec_xPro_hght = 200
tab2_dspec_yPro_wdth = tab2_dspec_wdth - tab2_dspec_xPro_wdth
tab2_dspec_yPro_hght = tab2_dspec_xPro_hght
tab2_aia_wdth = 400 - 5
tab2_aia_hght = 400 - 10
tab2_vla_wdth = 400 - 70
tab2_vla_hght = tab2_aia_hght
tab2_dspec_thumb_wdth = 400 - 70
tab2_dspec_thumb_hght = tab2_vla_hght

tab3_aia_submap_wdth = 720
tab3_aia_submap_hght = 720
tab3_dspec_small_wdth = 500
tab3_dspec_small_hght = 180
WebGL = False

config_plot = {'datadir': {'database': '/Users/fisher/Desktop/work/2016/NJIT/2014-11-01/test/database/',
                           'event_id': 'SUN01_20141101.163940-164700.50ms/',
                      'event_specfile': 'SUN01_20141101.163940-164700.50ms.cal.ms.spec.npz',
                      'fits_LOCL': 'Synthesis_Image/local/', 'fits_GLOB': 'Synthesis_Image/global/', 'J2000': 'J2000/',
                      'dspecDF': 'dspecDF-save', 'fits_GLOB_init': 'QLook/static/Synthesis_Image_fits_GLOB_init.fits',
                           'fits_LOCL_init': 'QLook/static/Synthesis_Image_fits_LOCL_init.fits'}, 'plot_config': {
    'tab_QLook': {'dspec_wdth': tab1_dspec_wdth, 'dspec_hght': tab1_dspec_hght, 'StrID_DataTb_BUT_wdth': 200,
             'dspec_Ctrl_widget_wdth': 200, 'StrID_DataTb_wdth': tab1_StrID_DataTb_wdth,
             'StrID_DataTb_hght': tab1_StrID_DataTb_hght},
    'tab_FSview_base': {'dspec_wdth': 1290, 'dspec_hght': 220, 'dspec_xPro_wdth': 900, 'dspec_xPro_hght': 200,
             'dspec_yPro_wdth': tab2_dspec_wdth - tab2_dspec_xPro_wdth, 'dspec_yPro_hght': tab2_dspec_xPro_hght,
             'aia_wdth': 400 - 5, 'aia_hght': 400 - 10, 'vla_wdth': 400 - 70, 'vla_hght': 400 - 10,
             'dspec_thumb_wdth': tab2_dspec_thumb_wdth, 'dspec_thumb_hght': tab2_dspec_thumb_hght},
    'tab_FSview_FitANLYS': {'aia_submap_wdth': 720, 'aia_submap_hght': 720, 'dspec_small_wdth': 500, 'dspec_small_hght': 180},
    'WebGL': False,
    'tab_ToClean': {'dspec_wdth': tab1_dspec_wdth, 'dspec_hght': tab1_dspec_hght - 100, 'dspec_xPro_wdth': 800,
                   'dspec_xPro_hght': 200, 'dspec_yPro_wdth': tab1_dspec_wdth - 800 - 150,
                   'dspec_yPro_hght': tab2_dspec_xPro_hght, 'aia_wdth': 400 - 5, 'aia_hght': 400 - 10,
                   'vla_wdth': 400 - 70, 'vla_hght': 400 - 10, 'dspec_thumb_wdth': tab2_dspec_thumb_wdth,
                   'dspec_thumb_hght': tab2_dspec_thumb_hght, 'input_tCLN_wdth': 500, 'input_tCLN_hght': 250,
                   'tab2_Div_tCLN_wdth': 800, 'tab2_Div_tCLN_hght': 300}}}

with open('config.json', 'w') as fp:
    json.dump(config_plot, fp)
