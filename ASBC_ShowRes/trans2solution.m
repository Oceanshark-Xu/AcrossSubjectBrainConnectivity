function cornew = trans2solution(corori,T1,T2)

mni = T1*[corori(:,1) corori(:,2) corori(:,3) ones(size(corori,1),1)]';
mni = mni';
mni(:,4) = [];

cornew = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T2))';
cornew(:,4) = [];



