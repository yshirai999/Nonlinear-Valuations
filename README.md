# Nonlinear Valuations

- In classical finance, valuations are discounted expectations of future cash flows. Equivalently, valuations are expectations of future cash flows under a risk neutral measure.

- In continuous time, the tower property of expectations allows one to define valuations as solutions of a certain PDE, such as the Black-Scholes equation.
  - When the value is allowed to jump we need to solve a partial integro differential equation (PIDE), where the additional integral term defines the expected jump of the price over the next instant

- For many practical purposes it is useful to have a range of possible prices instead of a single valuation operator. For instance:
  - To compute the market impact of a trade, one needs to develop models for bid and ask prices, rather than mid prices
  - To assess when is the optimal time to enter or exit the market, one needs a quantitative model for upper and lower values of the asset

- In this case, the PIDE for the upper and lower valuation bounds becomes a semilinear PIDE

- Appropriate models for the ranges of prices are thus of fundamental importance in finance, since:
  - A valuation range that is too large likely implies that we miss a chance to trade
  - If too small on the other hand we may not properly take into account all the possible scenarios.

- In this research, we
  - propose a few continuous time nonlinear valuation models within the framework of continuous time conic finace (aka: dynamic spectral risk measures);
  - calibrate the models to options on SPY
  - estimate the models for each of 9 sector ETFs and SPY based on their historical 5-day high and low prices;
  - Implement a trading strategy that maximizes the lower value of the portfolio.

For all the mathematical details, as well as the calibration and estimation part, see the accompanying paper [Continuous Time Conic Finance](https://www.aimsciences.org/article/doi/10.3934/fmf.2023021)

## Data

- All the calibration and estimation scripts are based on data that must be in the Data folder

- Daily d COB ata on Bilateral Gamma parameters is available for the period 01/02/2008 to 12/31/2020 and in the four struct files BGP1-BGP4.mat for several assets and indexes

- Data on options on SPY and ETFs can be downloaded from WRDS using the respective functions in the DataProcessing folder
  - These functions will prompt the user to input  WRDS username and password before the download begins

- For options on SPY data, the Data_Option.m function can be used
  - Data used in this project on option is limited to 1M maturity options, that are at most 30% OTM
  - The struct file created with the function Data_Option.m reports:
    - secid of the underlying asset (SPY)
    - Current date
    - Expiration date
    - Asset price
    - Strike price
    - Bid, Mid, Ask

- For ETFs data, the Data_VIX.m function should be used
  - Once ETFs data is downloaded, run the Data.m script to associate BG parameters to the ETF price for each day considered

## Calibration

- Calibration is performed based on options on SPY and for three different models:
  - BG2BG, which assumes that the upper, mid and lower prices follow exponential bilateral gamma (BG) processes
  - CGMY2CGMY, which assumes that the upper, mid and lower prices follow exponential CGMY processes
  - GcgammaDouble, which assumes that mid price is an exponential BG process and the upper and lower prices are exponential Levy processes with Levy measures: $\bar{\nu}(dy) = (1+\bar{\psi}_{\Gamma}(y))\nu(dy), \quad \underline{\nu}(dy) = (1+\underline{\psi}_{\Gamma}(y))\nu(dy) $
  where $\nu(dy)$ is the BG density, and $\underline{\psi}$ and $\overline{\psi}$ are given in proposition 4.1 of the accompanying paper

- The main goals of this calibration exercise are:
  - to assess if continuous time conic finance upper and lower prices can match relative bid-ask spreads across strikes; and
  - to infer, based on the different relative bid-ask spreads across strikes, which events the market is more uncertain about.

- Our results show that:
  - the GcgammaDouble model outperforms BG2BG and CGMY2CGMY, although the calibration is far slower since one cannot use the FFT method
  - when calibrating the upper price process to the ask price of calls and the bid price of puts, the market focuses on scenarios in which losses are lower than expected
  - when calibrating the lower price process to the bid price of calls and the ask price of puts, the market focuses on scenarios in which gains are larger than expected
  - For a more complete description of this results, see section 5.3 of the accompanying paper.

## Estimation
