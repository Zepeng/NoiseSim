int GenerateNoise(int seed,char* filename)
{
    std::cout << "start" << std::endl;
    TFile* wf_file = new TFile("digi.root");
    std::vector<double> *oversampling = 0;
    TTree* wftree = (TTree*)wf_file->Get("waveformTree");
    wftree->SetBranchAddress("OverSamplingAmp", &oversampling);
    wftree->GetEntry(2);
    //for(size_t iter = 0; iter < oversampling->size(); iter++)
    //    std::cout << oversampling->at(iter) << std::endl;
    
    TFile* noise_file = new TFile("noise_Aldo.root","read");

    TH1F* noise_hist = (TH1F*)noise_file->Get("h2");
    /*std::vector<float> noise_freq;
    std::vector<float> noise_mag;
    for(int i = 0; i < noise_hist->GetNbinsX(); i++)
    {
        if(noise_hist->GetBinCenter(i+1)>1800 || noise_hist->GetBinCenter(i+1)<200) continue;
        noise_freq.push_back(noise_hist->GetBinCenter(i+1));
        noise_mag.push_back(noise_hist->GetBinContent(i+1));
    }
    noise_file->Close();
    */
    TFile *f = TFile::Open(filename,"RECREATE");
    std::vector<float> noise_tdomain;
    std::vector<float> noise_filter;
    std::vector<float> noise_current;
    std::vector<float> noise_int;
    TTree *t = new TTree("noiselib","Noise Library");
    t->Branch("noise_int", &noise_current);
    //t->Branch("noise_int", &noise_int);
    //t->Branch("noise_filter", &noise_filter);
    //t->Branch("noise_tdomain", &noise_tdomain);

    int ntick = 220000;//1/fSamplingInterval*1.0e6;
	//Create Noise spectrum in frequency
	Double_t *noise_re = new Double_t[ntick];
	Double_t *noise_im = new Double_t[ntick];
	double pval = 0;
	double phase = 0;
	double rnd[2] = {0.};
    TRandom3* fRand = new TRandom3(seed);
	//width of frequency bin in kHz, add noise up to OverSamplingRatio*default sampling frequency 
	double binWidth = 25*1000./0.5/double(ntick);
	for(int lib_len = 0; lib_len < 100; lib_len++){
        for (int i = 0; i < ntick; i++)
	    {
		    double local_freqency = (i*binWidth + 0.5*binWidth)*1000;
		    //The noise spectrum stored in ROOT file has only frequency between ~30 kHZ to 1 MHz. Here, it's assumed white noise
		    //at frequencies higher than 2 MHz.
		    if(local_freqency < noise_hist->GetBinCenter(1))
            {
                pval = noise_hist->GetBinContent(1);
            }
            else if(local_freqency > noise_hist->GetBinCenter(110))
            {
                pval = noise_hist->GetBinContent(110);
            }
            else//Generate noise based on the noise spectrum in file
            {
                noise_hist->Interpolate(local_freqency);
            }
            //convert noise from nV to number of electrons using gain of 25 mV/fC
		    pval = pval/(25.*1000)*6240;
		    rnd[0] = fRand->Rndm();
		    rnd[1] = fRand->Rndm();
		    //Randomize noise magnitude by 10%
		    pval = pval*(0.9+0.2*rnd[0]);
		    //Random phase angle
		    phase = rnd[1]*2.*TMath::Pi();
		    if(i==0) pval = 0;//Turn off 0 frequency
		    noise_re[i] = pval*TMath::Cos(phase);
		    noise_im[i] = pval*TMath::Sin(phase);
	    }
	    //Obtain time spectrum from frequency spectrum
	    //noise_tdomain.resize(ntick,0.0);
        noise_tdomain.clear();
	    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &ntick, "C2R M K");
	    fftc2r->SetPointsComplex(noise_re, noise_im);
	    fftc2r->Transform();
	    double factor = 150000./sqrt((double)ntick);//1.8 added to normalize noise magnitude to have 200 e- noise per data point on current waveform of 2 MHz sampling rate. 
	
        for(int i = 0; i < 25*2000; i++)
        //for(int i=0; i < oversampling->size()/2; i++)
	    {
            //Save noise only waveform
            noise_tdomain.push_back(factor*fftc2r->GetPointReal(i+10000, false));
            //Save signal only waveform, the input waveform has a 100 MHz oversampling rate, need to dowmsample to 50 MHz.
            //noise_tdomain.push_back(oversampling->at(i*2));
            //save signal+noise waveform.
            //noise_tdomain.push_back(factor*fftc2r->GetPointReal(i+10000, false) + oversampling->at(i*2)) ;
	    }
        std::vector<double> wf_current;
        for(int iter = 0; iter < noise_tdomain.size() - 1; iter++)
            wf_current.push_back(noise_tdomain[iter+1] - noise_tdomain[iter]/2.);//Reduce the noise magnitude to 100 e- for sensitivity paper.
        double GAIN  = 4.206410398e+07;// 1.270891926e+09;
    
        static float xv[6] = {0}, yv[6] = {0};
        noise_filter.clear();
        for (int iter = 0; iter < wf_current.size();iter++)
        {
            xv[0] = xv[1]; xv[1] = xv[2]; xv[2] = xv[3]; xv[3] = xv[4]; xv[4] = xv[5];
            xv[5] = wf_current[iter]/ GAIN;
            yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5]; 
            yv[5] =   (xv[0] + xv[5]) + 5 * (xv[1] + xv[4]) + 10 * (xv[2] + xv[3])
            + (  0.7921866153 * yv[0]) + ( -4.1468660635 * yv[1])
            + (  8.6861252671 * yv[2]) + ( -9.1003509681 * yv[3])
            + (  4.7689043885 * yv[4]);
            /*yv[0] = yv[1]; yv[1] = yv[2]; yv[2] = yv[3]; yv[3] = yv[4]; yv[4] = yv[5];
            yv[5] =   (xv[0] + xv[5]) + 5 * (xv[1] + xv[4]) + 10 * (xv[2] + xv[3])
            + (  0.8900484892 * yv[0]) + ( -4.5543108339 * yv[1])
            + (  9.3224592095 * yv[2]) + ( -9.5421766696 * yv[3])
            + (  4.8839797796 * yv[4]);
            */
            noise_filter.push_back(yv[5]);
        }
        
        noise_current.clear();
        noise_int.clear();
        for(int noise_iter = 10; noise_iter < noise_filter.size()/25 - 10; noise_iter++)
            noise_current.push_back(noise_filter[noise_iter*25+70]*25);
        for(int noise_iter = 0; noise_iter < noise_current.size() - 1; noise_iter++)
            if(noise_iter == 0)
                noise_int.push_back(noise_current[noise_iter]);
            else
                noise_int.push_back(noise_current[noise_iter] + noise_int[noise_iter - 1]);
        t->Fill();
        f->Write();
    }
    return 0;
}
