function imgdef = translateImage(tx, ty, img)

A = [1, 0, 0; 0, 1, 0; tx, ty, 1];
tform = affine2d(A);

imgdef = imwarp(img, tform);
end