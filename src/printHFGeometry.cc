#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <iomanip>
#include <iostream>
#include <vector>



class PrintHFGeometry : public edm::one::EDAnalyzer<> {
public:
  explicit PrintHFGeometry(edm::ParameterSet const&)
      : caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()) {}

  void analyze(edm::Event const&, edm::EventSetup const& iSetup) override {
    auto const& caloGeometry = iSetup.getData(caloGeomToken_);

    auto const& detIds =
        caloGeometry.getValidDetIds(DetId::Hcal, HcalForward);

    std::cout << "Found " << detIds.size() << " HF cells\n";
    std::cout << std::fixed << std::setprecision(6);

    std::cout
        << "ieta iphi depth"
        << " eta phi"
        << " frontArea_cm2 backArea_cm2 avgArea_cm2"
        << " deta dphi"
        << "\n";

    for (auto const& id : detIds) {
      HcalDetId hid(id);

      auto const cell = caloGeometry.getGeometry(id);
      if (!cell.get()) {
        std::cout << "No CaloCellGeometry for rawId=" << id.rawId() << "\n";
        continue;
      }

      GlobalPoint const pos = cell->getPosition();

      double const frontArea = quadArea(cell->getCorners(), 0, 1, 2, 3);
      double const backArea  = quadArea(cell->getCorners(), 4, 5, 6, 7);
      double const avgArea   = 0.5 * (frontArea + backArea);

      std::cout
          << std::setw(4) << hid.ieta()
          << " " << std::setw(4) << hid.iphi()
          << " " << std::setw(2) << hid.depth()
          << " " << std::setw(12) << pos.eta()
          << " " << std::setw(12) << pos.phi()
          << " " << std::setw(14) << frontArea
          << " " << std::setw(14) << backArea
          << " " << std::setw(14) << avgArea
          << " " << std::setw(10) << cell->etaSpan()
          << " " << std::setw(10) << cell->phiSpan()
          << "\n";
    }
  }

private:
  using Corners = CaloCellGeometry::CornersVec;

  static GlobalVector diff(GlobalPoint const& a, GlobalPoint const& b) {
    return GlobalVector(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
  }

  static double triangleArea(GlobalPoint const& a,
                             GlobalPoint const& b,
                             GlobalPoint const& c) {
    GlobalVector const ab = diff(b, a);
    GlobalVector const ac = diff(c, a);
    return 0.5 * ab.cross(ac).mag();
  }

  static double quadArea(Corners const& c,
                         unsigned int i0,
                         unsigned int i1,
                         unsigned int i2,
                         unsigned int i3) {
    // Split quadrilateral into two triangles.
    return triangleArea(c[i0], c[i1], c[i2]) +
           triangleArea(c[i0], c[i2], c[i3]);
  }

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
};

DEFINE_FWK_MODULE(PrintHFGeometry);
