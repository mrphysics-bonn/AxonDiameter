% read in bvecs and bvals from grad file output from MRtrix
% should be able to use bval and bvec instead
grad=importdata('grad',' ',1);
bval=grad.data(:,4);
bvec=grad.data(:,1:3);

gd=niftiread('grad_dev.nii.gz');
dims=size(gd);
dims=dims(1:3);
gd=reshape(gd,[dims,3,3]);
gd=permute(gd,[5,4,1,2,3]);

bidxs{1}=bval==6000;
bidxs{2}=bval==30450;

meanbval=cell(2,1);
for bidx=1:length(bidxs)
    meanbval{bidx}=zeros(dims);
    for n=1:numel(meanbval{bidx})
        [i,j,k]=ind2sub(dims,n);
        L=eye(3)+squeeze(gd(:,:,i,j,k));

        bvecL=bvec(bidxs{bidx},:)*L;

        meanbval{bidx}(n)=mean(bval(bidxs{bidx}).*sum(bvecL.^2,2));
    end
end
