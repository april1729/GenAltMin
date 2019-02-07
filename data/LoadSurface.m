function LoadSurface

fname = 'data/al_1';

[pt trg] = ReadOFF([fname '.off']);

figure;
ViewMesh(pt,trg);