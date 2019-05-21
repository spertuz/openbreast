function showcase(dataset, id_no, study_no)

i = find((dataset.id==id_no)&(dataset.study==study_no));
for n = 1:length(i)
    figure
    info = getinfo(dataset.path{i(n)});
    im = ffdmRead(dataset.path{i(n)}, info);
    im = imresize(im, .1);    
    imshow(im)
    title(sprintf('%s%s',info.side, info.view))
end
    