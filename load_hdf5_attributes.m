function opt = load_hdf5_attributes(fdir, shotname)

dev_list = 'LHG';

opt.shotname = shotname;
shotnum = str2double(shotname);

for i = 1:3
    fname = sprintf('%s/ECEI.%06d.%cFS.h5', fdir, shotnum, dev_list(i));

    if exist(fname,'file') && (shotnum >= 7143 && shotnum < 9740)
        opt.Bt          = h5readatt(fname,'/ECEI','TFcurrent')*0.0995556;
        TriggerTime  = h5readatt(fname,'/ECEI','TriggerTime');
        opt.toffset(i)  = TriggerTime(1);
%         HarmonicMode    = cell2mat(h5readatt(fname,'/ECEI','HarmonicMode'));
        % opt.harmonic = str2double(HarmonicMode(2));
        opt.harmonic = [2,2,2];
        SampleRate       = h5readatt(fname,'/ECEI','SampleRate')*1000.0;
        opt.Fs(i) = SampleRate(1);
        opt.LO(i)       = h5readatt(fname,'/ECEI','LoFreq');
        opt.vzf(i)      = h5readatt(fname,'/ECEI','ZoomCodeV');
        opt.LensFocus(i) = h5readatt(fname,'/ECEI','LensFocus');
        opt.LensZoom(i) = h5readatt(fname,'/ECEI','LensZoom');
        % opt.voffset(i)  = 0;
    end

    if exist(fname,'file') && (shotnum > 9740 && shotnum < 12273)
        opt.Bt          = h5readatt(fname,'/ECEI','TFcurrent')*0.0995556;
        TriggerTime  = h5readatt(fname,'/ECEI','TriggerTime');
        opt.toffset(i)  = TriggerTime(1);
        opt.harmonic = [2,2,2];
        SampleRate       = h5readatt(fname,'/ECEI','SampleRate')*1000.0;
        opt.Fs(i) = SampleRate(1);
        opt.LO(i)       = h5readatt(fname,'/ECEI','LoFreq');
        opt.LensFocus(i) = h5readatt(fname,'/ECEI','LensFocus');
        opt.LensZoom(i) = h5readatt(fname,'/ECEI','LensZoom');
        % opt.vzf(i)      = h5readatt(fname,'/ECEI','ZoomCodeV');
        % opt.voffset(i)  = 0; 
    end

    if exist(fname,'file') && (shotnum > 12272 && shotnum < 30000)
        opt.Bt          = h5readatt(fname,'/ECEI','TFcurrent')*0.0995556;
        opt.TriggerTime{i}  = h5readatt(fname,'/ECEI','TriggerTime');
        opt.toffset(i)  = opt.TriggerTime{i}(1);
        SampleRate       = h5readatt(fname,'/ECEI','SampleRate')*1000.0;

        opt.Fs(i) = SampleRate(1);
        opt.LO(i)       = h5readatt(fname,'/ECEI','LoFreq');
        opt.LensFocus(i) = h5readatt(fname,'/ECEI','LensFocus');
        opt.LensZoom(i) = h5readatt(fname,'/ECEI','LensZoom');
        Mode = char(h5readatt(fname,'/ECEI','Mode'));
        if Mode == 'O'
            opt.harmonic(i) = 1;     
        elseif Mode == 'X'
            opt.harmonic(i) = 2;     
        end
    end
end

end