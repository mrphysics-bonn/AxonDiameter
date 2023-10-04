# Imports

import ants
import antspynet
import argparse
import nibabel as nib

def main(args):

    nifti = nib.load(args.in_file)
    ants_img = ants.image_read(args.in_file)

    mask = antspynet.brain_extraction(ants_img, modality="t1")
    denoise = ants.denoise_image(ants_img,mask=mask)
    n4 = ants.n4_bias_field_correction(denoise,mask=mask)

    thresh = 0.2
    mask = mask.numpy()
    mask[mask>thresh] = 1
    mask[mask<=thresh] = 0
    t1_preprocessed = n4.numpy()*mask

    nifti_img = nib.Nifti1Image(t1_preprocessed, affine=nifti.affine, header=nifti.header)
    nib.save(nifti_img, args.out_file+'.nii.gz')

    if args.mask:
        nifti_img = nib.Nifti1Image(mask, affine=nifti.affine, header=nifti.header)
        nib.save(nifti_img, args.out_file+'_mask.nii.gz')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='T1 image preprocessing (brain extraction denoising and n4 bias correction)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_file', type=str, help='Input T1 file')
    parser.add_argument('out_file', type=str, help='Output filename for preprocessed T1 data (without .nii.gz ending).')
    parser.add_argument('-m', '--mask', action='store_true', help='Output the brain mask.')

    args = parser.parse_args()

    main(args)
