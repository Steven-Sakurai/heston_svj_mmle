#include <iostream>
#include <iomanip>
#include "loglike.hpp"

int main() {
	cout << setprecision(16) << setw(16);
	cout << "test: " << endl;
	
	stdVec par = {0.05, 2, 0.2, 0.25, -0.8};
	stdVec y = {4.6052,4.6195,4.5787,4.5915,4.6029,4.5893,4.582,4.6064,4.6063,4.591,4.603,4.6037,4.6085,4.6084,4.6098,4.6184,4.5885,4.5867,4.5745,4.5805,4.5992,4.6297,4.6185,4.6292,4.6224,4.614,4.5933,4.5941,4.6392,4.6403,4.6247,4.6329,4.6567,4.6437,4.6464,4.6546,4.642,4.6428,4.6205,4.6622,4.6615,4.6971,4.6814,4.6613,4.628,4.6256,4.5843,4.6027,4.629,4.618,4.6309,4.6188,4.6193,4.6032,4.6018,4.6181,4.6088,4.5947,4.5965,4.6116,4.575,4.5864,4.596,4.6058,4.6073,4.6021,4.5985,4.6157,4.5865,4.585,4.5715,4.5654,4.5689,4.579,4.5927,4.5798,4.5738,4.5562,4.5999,4.6129,4.6286,4.631,4.6352,4.6328,4.6573,4.6416,4.6114,4.6169,4.6244,4.6429,4.6373,4.6484,4.6766,4.6824,4.6743,4.6623,4.6725,4.6676,4.6729,4.7033,4.6877,4.707,4.7061,4.6905,4.6757,4.6634,4.6915,4.6806,4.6799,4.6925,4.6701,4.6605,4.6784,4.6845,4.6828,4.687,4.6537,4.6304,4.6326,4.6548,4.6362,4.6061,4.5576,4.5512,4.5874,4.579,4.5851,4.5941,4.5829,4.6218,4.5748,4.6018,4.5787,4.6437,4.6679,4.6652,4.6579,4.6825,4.71,4.6888,4.6866,4.6974,4.6899,4.6904,4.6727,4.6804,4.7244,4.7448,4.7542,4.7305,4.719,4.7526,4.7357,4.7559,4.7678,4.7793,4.8134,4.8067,4.8069,4.7959,4.8067,4.8084,4.8033,4.8021,4.8188,4.8122,4.7911,4.7995,4.7996,4.8115,4.821,4.8326,4.8344,4.8381,4.8476,4.8521,4.8578,4.8689,4.8561,4.858,4.8869,4.8646,4.8563,4.8511,4.848,4.8451,4.7824,4.7788,4.7882,4.7941,4.7752,4.7719,4.7632,4.7681,4.7586,4.7642,4.7619,4.7417,4.7364,4.7443,4.7217,4.7325,4.7334,4.7443,4.7466,4.7273,4.7297,4.7169,4.7146,4.6956,4.6872,4.6644,4.6561,4.6162,4.6529,4.6778,4.6656,4.6523,4.6221,4.5825,4.5816,4.5715,4.5725,4.569,4.569,4.5509,4.5684,4.5769,4.5636,4.5926,4.5784,4.5924,4.5905,4.5923,4.595,4.5886,4.609,4.6287,4.6396,4.618,4.6162,4.6454,4.6605,4.6546,4.6242,4.6372,4.5972,4.6103,4.615,4.6331,4.6293,4.6427,4.6301,4.6334,4.6068,4.6178,4.592,4.6094,4.6133,4.5997,4.5712,4.5694,4.5852,4.5599,4.5884,4.5759,4.5904,4.5801,4.5766,4.58,4.5638,4.5749,4.5569,4.5987,4.6015,4.5897,4.5879,4.5962,4.5656,4.5217,4.5274,4.464,4.4698,4.4618,4.4612,4.4584,4.4537,4.4815,4.4822,4.4894,4.4692,4.4862,4.5311,4.5289,4.547,4.577,4.5713,4.5831,4.6083,4.5963,4.6083,4.598,4.5954,4.6259,4.6034,4.5733,4.534,4.5589,4.6056,4.6041,4.5904,4.5718,4.538,4.563,4.5473,4.5257,4.5386,4.5555,4.5909,4.6094,4.6012,4.5927,4.6149,4.6382,4.6327,4.6238,4.555,4.5432,4.517,4.5342,4.4924,4.5109,4.5181,4.5316,4.5378,4.5874,4.5963,4.6473,4.6603,4.633,4.6252,4.6112,4.5841,4.5831,4.5835,4.5839,4.5929,4.5882,4.5792,4.5682,4.5664,4.6032,4.6207,4.6652,4.6467,4.6373,4.665,4.7005,4.6906,4.6998,4.7145,4.7438,4.7753,4.7574,4.7439,4.738,4.725,4.7428,4.7052,4.7342,4.7535,4.754,4.7457,4.744,4.7037,4.7059,4.7032,4.7175,4.7537,4.7603,4.7913,4.7946,4.7959,4.8241,4.8097,4.8231,4.8453,4.8649,4.8781,4.8756,4.8826,4.889,4.874,4.8618,4.8768,4.8895,4.9043,4.9138,4.9148,4.9059,4.9101,4.8826,4.898,4.9217,4.9383,4.925,4.9288,4.9463,4.9336,4.9777,4.9858,4.9795,4.9876,4.9852,4.987,4.9716,4.9756,4.9277,4.8948,4.9214,4.9127,4.9193,4.9276,4.9246,4.9378,4.9315,4.926,4.9419,4.9441,4.9512,4.953,4.9466,4.9632,4.9777,4.9422,4.9667,4.9845,4.9638,4.9667,4.9771,4.981,4.9865,4.9908,4.9857,4.9615,4.9576,4.9365,4.9226,4.9001,4.9119,4.8879,4.8929,4.9116,4.9127,4.9283,4.9191,4.923,4.9187,4.9041,4.8873,4.8825,4.9179,4.8968,4.8656,4.8793,4.8643,4.9059,4.8743,4.8903,4.8837,4.865,4.856,4.8368,4.8142,4.7996,4.8112,4.8056,4.8043,4.8174,4.8015,4.7893,4.7978,4.7966,4.788,4.7843,4.7961,4.7913,4.811,4.8157,4.8109,4.7987,4.8029,4.7964,4.7974,4.7823,4.7504,4.7606,4.7514,4.742,4.7646,4.7872,4.7983,4.7822,4.7756,4.7745,4.765,4.7954,4.8004,4.7827,4.7217,4.7302,4.6964,4.7101,4.6951,4.7207,4.7079,4.7523,4.7407,4.764,4.7371,4.7372,4.7373,4.7301,4.7146,4.704,4.6946,4.6986,4.709,4.6845,4.6979,4.7129,4.6833,4.6687,4.6598,4.6065,4.5957,4.6153,4.6523,4.6337,4.6034,4.5969,4.5874,4.5677,4.5452,4.5275,4.5402,4.5158,4.5428,4.5638,4.582,4.6063,4.6138,4.6207,4.5997,4.6011,4.5799,4.5716,4.5763,4.6399,4.6547,4.6411,4.6254,4.6148,4.5792,4.5461,4.576,4.5368,4.5129,4.5239,4.5179,4.4969,4.5124,4.5325,4.5702,4.5341,4.5194,4.4989,4.5311,4.5328,4.5574,4.5493,4.5501,4.5667,4.5651,4.5824,4.57,4.5457,4.5416,4.5574,4.5689,4.5762,4.5676,4.5661,4.5729,4.5546,4.5419,4.5541,4.5635,4.5405,4.5826,4.5835,4.6113,4.5903,4.6011,4.6047,4.6187,4.6185,4.6193,4.629,4.6652,4.671,4.7009,4.7096,4.6916,4.7139,4.7414,4.746,4.7522,4.7811,4.7708,4.7703,4.7539,4.7411,4.7528,4.7688,4.7316,4.7121,4.6964,4.7095,4.6898,4.6791,4.6745,4.6705,4.6667,4.6721,4.6786,4.6886,4.6864,4.6897,4.714,4.7519,4.7831,4.7575,4.7635,4.7309,4.7541,4.7594,4.7791,4.7595,4.7548,4.7594,4.7651,4.7299,4.7683,4.7447,4.695,4.7208,4.6983,4.7071,4.6956,4.7071,4.7265,4.7176,4.7091,4.6937,4.6892,4.6802,4.6678,4.6558,4.6453,4.636,4.623,4.6015,4.6304,4.6118,4.6479,4.6394,4.6307,4.6257,4.6392,4.6183,4.6259,4.6168,4.6409,4.6593,4.618,4.6154,4.5619,4.5584,4.4888,4.4723,4.4753,4.4617,4.5047,4.5055,4.5094,4.536,4.5299,4.5258,4.5097,4.4978,4.49,4.4943,4.4963,4.5078,4.4941,4.4693,4.459,4.4485,4.4423,4.4241,4.4559,4.4565,4.476,4.4767,4.5098,4.4955,4.5311,4.509,4.5379,4.5181,4.5231,4.554,4.5585,4.574,4.5735,4.5534,4.5473,4.5274,4.5229,4.5361,4.5406,4.5764,4.5583,4.5699,4.5641,4.5667,4.594,4.6071,4.6248,4.639,4.6648,4.6459,4.6458,4.6461,4.6268,4.6308,4.6607,4.6609,4.6617,4.6698,4.652,4.6462,4.6273,4.6038,4.6123,4.6142,4.5955,4.6163,4.6572,4.6232,4.627,4.6342,4.6314,4.6356,4.6035,4.591,4.6158,4.6172,4.5999,4.5724,4.6006,4.5914,4.5387,4.5547,4.5247,4.5173,4.5024,4.4981,4.4738,4.4999,4.4847,4.4483,4.4604,4.4556,4.4762,4.4777,4.5146,4.5307,4.5614,4.5596,4.5545,4.5447,4.5465,4.5442,4.5437,4.5634,4.5768,4.5867,4.5623,4.5408,4.5724,4.5846,4.5649,4.5763,4.5475,4.5435,4.5611,4.555,4.5247,4.5016,4.506,4.5094,4.5096,4.487,4.5138,4.4824,4.4475,4.4644,4.5183,4.5248,4.546,4.5685,4.5679,4.5403,4.5523,4.5386,4.5225,4.5179,4.5253,4.5393,4.5323,4.5108,4.487,4.4972,4.5224,4.5202,4.5413,4.517,4.5246,4.5249,4.5252,4.5174,4.5129,4.48,4.4917,4.4644,4.4693,4.4444,4.4317,4.4375,4.4209,4.4234,4.4193,4.3865,4.4081,4.4153,4.4679,4.4551,4.4624,4.4571,4.4852,4.4723,4.4684,4.5007,4.4948,4.4951,4.4647,4.5027,4.5297,4.5196,4.538,4.5155,4.5318,4.5441,4.5544,4.5344,4.5448,4.5347,4.519,4.525,4.5317,4.5261,4.5447,4.5101,4.5167,4.5058,4.4963,4.4995,4.4983,4.5103,4.52,4.526,4.5209,4.517,4.494,4.5036,4.4785,4.4944,4.4831,4.5036,4.5216,4.4673,4.4529,4.4557,4.45,4.3997,4.3951,4.4295,4.4136,4.4128,4.4286,4.4483,4.4556,4.4525,4.4284,4.4361,4.418,4.411,4.3878,4.4096,4.4256,4.4389,4.4292,4.46,4.4283,4.4233,4.4368,4.4432,4.4617,4.4659,4.4488,4.427,4.4537,4.4828,4.482,4.4478,4.4143,4.4141,4.4007,4.3697,4.4301,4.4528,4.4661,4.4152,4.457,4.4693,4.4321,4.4622,4.4297,4.4201,4.4417,4.471,4.4649,4.4687,4.4927,4.4791,4.4992,4.5059,4.5354,4.5625,4.5609,4.5408,4.5594,4.5576,4.5785,4.575,4.577,4.5881,4.6154,4.6004,4.6023,4.5969,4.621,4.5955,4.6093,4.5741,4.5661,4.5774,4.5796,4.5806,4.6061,4.6108,4.6221,4.6142,4.6117,4.6378,4.6473,4.6517,4.6465,4.6648,4.6501,4.6391,4.6456,4.6575,4.6465,4.6336,4.6297,4.6499,4.6305,4.64};
    cout << loglike(par, 1000) << endl;
}
