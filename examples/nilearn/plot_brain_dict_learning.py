import os
import glob
from nilearn.image.image import check_niimg
from nilearn.masking import unmask
from nilearn_prox_operators import ProximalOperator
from fmri_dict_learning import ProximalfMRIMiniBatchDictionaryLearning
import sys
sys.path.append(os.path.join(os.environ["HOME"],
                             "CODE/FORKED/elvis_dohmatob/rsfmri2tfmri"))
from config import hcp_distro, data_dir, root

random_state = 42
n_components = 40
bcd_n_iter = 1
n_epochs = 1
dict_alpha = 10.
batch_size = 30
dataset = "IBC bold"

if dataset == "HCP":
    from datasets import fetch_hcp_task
    mask_img = os.path.join(data_dir, hcp_distro, "mask_img.nii.gz")
    zmaps = fetch_hcp_task(os.path.join(data_dir, hcp_distro))
    X = zmaps[zmaps.contrast_name == "STORY-MATH"].groupby( 
       "subject_id")["zmap"].apply(sum)
elif dataset == "IBC zmaps":
    X = sorted(glob.glob(os.path.join(root,
                                      "storage/store/data/ibc",
                                      "derivatives/sub-*/ses-*",
                                      "res_stats_hcp_*_ffx",
                                      "stat_maps/*.nii.gz")))
    mask_img = os.path.join("/home/elvis/mnt/32-bit-system/home/elvis/drago",
    "storage/store/data/ibc/derivatives/group/gm_mask.nii.gz")
elif dataset == "IBC bold":
    from datasets import load_imgs
    mask_img = os.path.join("/home/elvis/mnt/32-bit-system/home/elvis/drago",
    "storage/store/data/ibc/derivatives/group/gm_mask.nii.gz")
    X = glob.glob(os.path.join(root,
                               "storage/store/data/ibc/derivatives",
                               "masked_IBC/sub-*/ses-05/func",
                               "wrdcsub-*_bold.npy"))
    assert 0
    X = [img for Xs in load_imgs(X) for img in Xs]
else:
    raise NotImplementedError(dataset)

batch_size = min(batch_size, len(X))
mask_img = check_niimg(mask_img)
n_voxels = mask_img.get_data().sum()

prox = ProximalOperator(which="graph-net", affine=mask_img.affine, fwhm=2,
                        mask=mask_img.get_data().astype(bool), l1_ratio=.1,
                        kernel="gaussian", radius=1., positive=True)
model = ProximalfMRIMiniBatchDictionaryLearning(
    n_components=n_components, random_state=random_state, fit_algorithm="ridge",
    dict_penalty_model=-2, mask=mask_img, n_epochs=1, batch_size=batch_size,
    dict_alpha=dict_alpha, verbose=2)
model.fit(X)

import matplotlib.pyplot as plt
from nilearn.plotting import plot_stat_map

for c in [3, 7, 17, 37, 38]:
    if c == 3:
        continue
    plot_stat_map(unmask(model.components_[c], mask_img), display_mode="yz",
                  colorbar=False, black_bg=True)
    out_file = os.path.join("/home/elvis/CODE/FORKED/coordescendant/paper",
                            "figs/%s_comp%02i.png" % (dataset, c))
    plt.savefig(out_file, dpi=200, bbox_inches="tight", pad_inches=0)
    os.system("mogrify -trim %s" % out_file)
    print(out_file)
plt.show()
