from pygimli.physics.SIP import SIPSpectrum

sip = SIPSpectrum('sipexample.txt')
sip.showData()  # show amplitude and phase spectra
sip.showDataKK()  # show data along with Kramers-Kronig calculated (data check)
sip.showData(znorm=True)  # show normalized imaginary/real parts instead
if 0:  # Pelton approach: fit Cole-Cole model along with EM Cole-Cole term
    sip.fitCCEM()
else:  # determine permittivity to correct values with
    sip.removeEpsilonEffect()
    sip.fitColeCole(useCond=False)
# %%
sip.fitDebyeModel(showFit=True)
# %% create titles and plot data, fit and model
sip.showAll(save=True)
