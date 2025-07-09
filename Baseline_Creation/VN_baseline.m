function VN_baseline

[stopperCount,countMatrix,segments]  = VN_sample_baseline(0,0,0,200);

VN_animate_baseline([segments.starts; segments.ends]);

% scaled version, doubled the position coordninates of the vertices
csvdata         = [segments.ends];

% name the csv file
csvname         = ['VN_baseline' ];

save('VN_segments','segments')

% write csv file
dlmwrite(csvname,csvdata);

end