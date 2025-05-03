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

- The purpose of this research is threefold: we
  - propose a conitnuous time nonlinear valuation model;
  - calibrate it to options on SPY
  - estimate it for each of 9 sector ETFs and SPY based on their historical 5-day high and low prices;
  - Implement a trading strategy that maximizes the lower value of the portfolio.

For more mathematical details, as well as the calibration and estimation part, see the accompanying paper [Continuous Time Conic Finance](https://www.aimsciences.org/article/doi/10.3934/fmf.2023021)
