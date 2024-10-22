function createCPPHeader(grad, rewind, slewRate0, slewRate1, ...
                         NumArms, ...
                         FOV0, FOV1, FOV2, baseRes, Resolution, ...
                         dwellTime, nstep, maxAmp, filePath)
    fID = fopen(filePath, 'w');

    fprintf(fID,'#pragma once\n\n\n');


    fprintf(fID, '/* NumArms = %d\n', NumArms);
    fprintf(fID, '   Max gradient amplitude = %.2f mT/m\n', maxAmp.*10);
    fprintf(fID, '   Read out max slew rate = %.2f mT/m/ms\n', slewRate0/100);
    fprintf(fID, '   Rewind max slew rate = %.2f mT/m/ms\n', slewRate1/100);
    fprintf(fID, ['   FOV = %d cm;\n' ...
                  '   FOVCoef1 = %d;\n' ...
                  '   FOVCoef2 = %d\n'], ...
            FOV0, FOV1, FOV2);
    fprintf(fID, '   Base resolution = %d\n', baseRes);
    fprintf(fID, '   Inplane Resolution = %.2f mm */\n\n\n', Resolution);

    
    fprintf(fID, 'long lDwellTime = %d; // ADC dwell time, micro seconds\n', ...
            uint64(1e6*dwellTime));
    fprintf(fID, 'long lDownsampleRate = %d; // micro seconds\n\n', ...
            uint64(nstep));

    fprintf(fID, 'double dMaxGrad = %f;\n\n', maxAmp.*10);

    fprintf(fID, 'double X[%d] = {\n', size(grad, 1));
    for i=1:length(grad)-1
        fprintf(fID, '    %.16f,\n', grad(i, 1)./maxAmp);
    end
    fprintf(fID, '    %.16f\n', grad(end, 1)./maxAmp);
    fprintf(fID, '};\n\n');

    fprintf(fID, 'double Y[%d] = {\n', size(grad, 1));
    for i=1:length(grad)-1
        fprintf(fID, '    %.16f,\n', grad(i, 2)./maxAmp);
    end
    fprintf(fID, '    %.16f\n', grad(end, 2)./maxAmp);
    fprintf(fID, '};\n\n');

    fprintf(fID, 'double XRW[%d] = {\n', size(rewind, 1));
    for i=1:length(rewind)-1
        fprintf(fID, '    %.16f,\n', rewind(i, 1)./maxAmp);
    end
    fprintf(fID, '    %.16f\n', rewind(end, 1)./maxAmp);
    fprintf(fID, '};\n\n');

    fprintf(fID, 'double YRW[%d] = {\n', size(rewind, 1));
    for i=1:length(rewind)-1
        fprintf(fID, '    %.16f,\n', rewind(i, 2)./maxAmp);
    end
    fprintf(fID, '    %.16f\n', rewind(end, 2)./maxAmp);
    fprintf(fID, '};\n');


    fclose(fID);


    % Trajectory files for reconstruction
    fID = fopen('Results/GX.txt', 'w');
    
    for i=1:length(grad)
        fprintf(fID, '%.16f,\n', grad(i, 1)./maxAmp);
    end

    fclose(fID);

    fID = fopen('Results/GY.txt', 'w');
    
    for i=1:length(grad)
        fprintf(fID, '%.16f,\n', grad(i, 2)./maxAmp);
    end

    fclose(fID);
end