#include "GPD3D_PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"   // twopi, pi
#include "G4GenericMessenger.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"             // G4UniformRand()

#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>

// ---------------------------------------------------------
// RNG for Gaussian beam spot
// ---------------------------------------------------------
static std::default_random_engine generator;
static std::normal_distribution<double> gauss01(0.0, 1.0);

// ---------------------------------------------------------
// Helpers: cosine-law hemisphere around a normal (Lambertian)
// ---------------------------------------------------------
static void OrthonormalBasis(const G4ThreeVector& n, G4ThreeVector& e1, G4ThreeVector& e2)
{
  // n assumed unit
  const G4ThreeVector a = (std::fabs(n.z()) < 0.999) ? G4ThreeVector(0,0,1) : G4ThreeVector(0,1,0);
  e1 = (n.cross(a)).unit();
  e2 = (n.cross(e1)).unit();
}

static G4ThreeVector SampleCosineHemisphere(const G4ThreeVector& n_in)
{
  // Cosine-law on hemisphere: cos(theta)=sqrt(u)
  const G4double u1 = G4UniformRand();
  const G4double u2 = G4UniformRand();

  const G4double cosA = std::sqrt(u1);
  const G4double sinA = std::sqrt(1.0 - u1);
  const G4double beta = twopi * u2;

  G4ThreeVector e1, e2;
  OrthonormalBasis(n_in, e1, e2);

  const G4ThreeVector dir =
      (sinA * std::cos(beta)) * e1 +
      (sinA * std::sin(beta)) * e2 +
      (cosA) * n_in;

  return dir.unit();
}

// ---------------------------------------------------------
// Constructor / destructor
// ---------------------------------------------------------
GPD3D_PrimaryGeneratorAction::GPD3D_PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fGPS(new G4GeneralParticleSource()),
  // Sphere 
  m_useSphere(true),
  m_sphereRadius(130.*mm),
  m_sphereCenter(0.,0.,0.),
  m_posThetaMin(0.0), // radians (0 = forward along +Z; pi = backward along -Z)
  m_posThetaMax(CLHEP::pi), // radians (0 = forward along +Z; pi = backward along -Z)
  // Spectrum energy from CSV file
  m_useSpectrum(false), //default is to use campana .mac
  m_spectrumFile(""),
  // Beam controls (used when useSphere == false)
  m_beamSigma(2.*mm),
  m_beamCenter(0.,0.,0.),
  m_beamZ(100.*mm),
 // GPS macro mode: if true, do not override GPS settings in C++; rely on /gps commands from macro
  m_useGpsMacro(false),
  // UI
  m_genMessenger(nullptr)
{
  // Seed RNG once
  const auto seed = static_cast<unsigned>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()
  );
  generator.seed(seed);

  // Default GPS config (kept from your original)
  auto* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGPS->SetParticleDefinition(proton);
  fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(100 * GeV);
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0., -100., 10.) * mm);

  // Polarization fixed at 30 degrees in XY plane (kept as-is)
  const G4double angle = 30.0 * CLHEP::pi / 180.0;
  fGPS->GetCurrentSource()->SetParticlePolarization(G4ThreeVector(std::cos(angle), std::sin(angle), 0.0));

  // UI commands
  m_genMessenger = new G4GenericMessenger(this, "/gpd3d/gen/", "Primary generation controls");

  // Sphere generation (4pi isotropic environment)
  m_genMessenger->DeclareProperty("useSphere", m_useSphere,
    "If true: generate primaries on a sphere with cosine-law inward directions.");
  m_genMessenger->DeclarePropertyWithUnit("sphereRadius", "mm", m_sphereRadius,
    "Sphere radius (must fit inside world).");
  m_genMessenger->DeclarePropertyWithUnit("sphereCenter", "mm", m_sphereCenter,
    "Sphere center (x y z mm).");
  m_genMessenger->DeclareProperty("useGpsMacro", m_useGpsMacro,
    "If true: do not override GPS settings in C++; rely on /gps commands from macro.");

  // Beam controls (when useSphere == false)
  m_genMessenger->DeclarePropertyWithUnit("beamSigma", "mm", m_beamSigma,
    "Gaussian beam sigma (beam mode only).");
  m_genMessenger->DeclarePropertyWithUnit("beamZ", "mm", m_beamZ,
    "Beam source Z position (beam mode only).");
  m_genMessenger->DeclarePropertyWithUnit("beamCenter", "mm", m_beamCenter,
    "Beam center (x y z mm). For beam mode, x/y are used; z is ignored (use beamZ).");
  m_genMessenger->DeclarePropertyWithUnit("posThetaMin", "rad", m_posThetaMin,
    "Minimum polar angle (theta) for particle direction generation.");
  m_genMessenger->DeclarePropertyWithUnit("posThetaMax", "rad", m_posThetaMax,
    "Maximum polar angle (theta) for particle direction generation.");
  // Energy spectrum controls
  m_genMessenger->DeclareProperty("useSpectrum", m_useSpectrum,
    "If true: sample primary energy from a CSV spectrum (log-log interpolation between points)."
    " If false: use the mono-energy already configured in the GPS.");

  // Use a method instead of DeclareProperty so we can reload immediately when changed from macro.
  m_genMessenger->DeclareMethod("spectrumFile", &GPD3D_PrimaryGeneratorAction::SetSpectrumFile,
    "CSV file path with 2 columns: energy_MeV, weight (arbitrary units). Lines starting with # are ignored.");

}

GPD3D_PrimaryGeneratorAction::~GPD3D_PrimaryGeneratorAction()
{
  delete m_genMessenger;
  delete fGPS;
}

// ---------------------------------------------------------
// GeneratePrimaries
// ---------------------------------------------------------
void GPD3D_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  // Pure GPS-from-macro mode: do not override pos/dir/energy
  if (m_useGpsMacro) {
    fGPS->GeneratePrimaryVertex(anEvent);
    return;
  }

  // Optional: sample energy from spectrum every event
  if (m_useSpectrum && !m_specCum.empty()) {
    const G4double e = SampleSpectrumEnergy();
    fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(e);
  }

  // --- Sphere mode: 4pi isotropic environment injection ---
  if (m_useSphere) {
    // --- Uniform point on sphere, restricted to theta in [m_posThetaMin, m_posThetaMax] ---
    G4double tmin = m_posThetaMin;
    G4double tmax = m_posThetaMax;

    // Clamp and protect against bad user input
    if (tmin < 0.0) tmin = 0.0;
    if (tmax > CLHEP::pi) tmax = CLHEP::pi;
    if (tmax <= tmin) { tmin = 0.0; tmax = CLHEP::pi; } // fallback: full sphere

    // Sample cos(theta) uniformly between cos(tmax) and cos(tmin)
    const G4double u = G4UniformRand();
    const G4double v = G4UniformRand();

    const G4double cosTmin = std::cos(tmax); // note: theta increases => cos decreases
    const G4double cosTmax = std::cos(tmin);

    const G4double cosT = cosTmin + u * (cosTmax - cosTmin);
    const G4double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
    const G4double phi  = twopi * v;

    const G4ThreeVector rhat(sinT*std::cos(phi), sinT*std::sin(phi), cosT);
    const G4ThreeVector pos = m_sphereCenter + m_sphereRadius * rhat;

    // Inward normal toward center
    const G4ThreeVector n_in = (m_sphereCenter - pos).unit();

    // Cosine-law direction inward (correct for isotropic intensity crossing a surface)
    const G4ThreeVector dir = SampleCosineHemisphere(n_in);

    fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
    fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(pos);
    fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(dir);

    fGPS->GeneratePrimaryVertex(anEvent);
    return;
  }

  // --- Beam mode: Gaussian spot at fixed Z ---
  const G4double x = m_beamCenter.x() + gauss01(generator) * m_beamSigma;
  const G4double y = m_beamCenter.y() + gauss01(generator) * m_beamSigma;

  const G4ThreeVector pos(x, y, m_beamZ);

  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(pos);

  fGPS->GeneratePrimaryVertex(anEvent);
}

// ---------------------------------------------------------
// Spectrum handling
// ---------------------------------------------------------

void GPD3D_PrimaryGeneratorAction::SetSpectrumFile(const G4String& path)
{
  m_spectrumFile = path;
  LoadSpectrumFromCsv(m_spectrumFile);
}

void GPD3D_PrimaryGeneratorAction::LoadSpectrumFromCsv(const G4String& path)
{
  m_specE.clear();
  m_specW.clear();
  m_specSegA.clear();
  m_specCum.clear();

  if (path.empty()) return;

  std::ifstream in(path);
  if (!in) {
    G4Exception("GPD3D_PrimaryGeneratorAction::LoadSpectrumFromCsv",
                "GPD3D_SPECTRUM_0001",
                JustWarning,
                ("Cannot open spectrum CSV: " + path).c_str());
    return;
  }

  std::string line;
  while (std::getline(in, line)) {
    // trim
    const auto trim = [](const std::string& s) -> std::string {
      const auto b = s.find_first_not_of(" \t\r\n");
      if (b == std::string::npos) return "";
      const auto e = s.find_last_not_of(" \t\r\n");
      return s.substr(b, e - b + 1);
    };

    const auto startsWithHash = [](const std::string& s) -> bool {
      for (char ch : s) {
        if (ch == ' ' || ch == '\t') continue;
        return ch == '#';
      }
      return false;
    };

    line = trim(line);
    if (line.empty() || startsWithHash(line)) continue;

    // Accept comma or whitespace separated; also accept ';' as comma
    for (char& c : line) if (c == ';') c = ',';

    double e = 0.0;
    double w = 0.0;
    {
      std::stringstream ss(line);
      std::string a, b;
      if (!std::getline(ss, a, ',')) continue;
      if (!std::getline(ss, b, ',')) {
        // maybe whitespace separated
        ss.clear();
        ss.str(line);
        if (!(ss >> e >> w)) continue;
      } else {
        try {
          e = std::stod(trim(a));
          w = std::stod(trim(b));
        } catch (...) {
          continue;
        }
      }
    }

    if (e <= 0.0 || w <= 0.0) continue;
    m_specE.push_back(e * MeV); // CSV energy is interpreted as MeV
    m_specW.push_back(w);
  }

  if (m_specE.size() < 2) {
    G4Exception("GPD3D_PrimaryGeneratorAction::LoadSpectrumFromCsv",
                "GPD3D_SPECTRUM_0002",
                JustWarning,
               "Spectrum CSV must contain at least two (E_MeV, weight) points with positive values.");
    m_specE.clear();
    m_specW.clear();
    return;
  }

  // Sort by energy
  std::vector<std::size_t> idx(m_specE.size());
  for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(),
            [&](std::size_t a, std::size_t b){ return m_specE[a] < m_specE[b]; });

  std::vector<G4double> e2; e2.reserve(idx.size());
  std::vector<G4double> w2; w2.reserve(idx.size());
  for (auto i : idx) { e2.push_back(m_specE[i]); w2.push_back(m_specW[i]); }
  m_specE.swap(e2);
  m_specW.swap(w2);

  // Build per-segment power-law exponent "a" for log-log interpolation and cumulative integral.
  // Between (E0,W0) and (E1,W1): W(E) = k * E^a with a = ln(W1/W0)/ln(E1/E0)
  // Segment integral: ∫ W(E) dE from E0..E1
  const std::size_t nseg = m_specE.size() - 1;
  m_specSegA.resize(nseg);
  m_specCum.resize(nseg + 1);
  m_specCum[0] = 0.0;

  for (std::size_t i = 0; i < nseg; ++i) {
    const G4double E0 = m_specE[i];
    const G4double E1 = m_specE[i+1];
    const G4double W0 = m_specW[i];
    const G4double W1 = m_specW[i+1];

    if (E1 <= E0 || W0 <= 0.0 || W1 <= 0.0) {
      m_specSegA[i] = 0.0;
      continue;
    }

    const G4double a = std::log(W1/W0) / std::log(E1/E0);
    m_specSegA[i] = a;

    // k = W0 / E0^a
    const G4double k = W0 / std::pow(E0, a);
    G4double segInt = 0.0;
    if (std::fabs(a + 1.0) > 1e-14) {
      segInt = k / (a + 1.0) * (std::pow(E1, a + 1.0) - std::pow(E0, a + 1.0));
    } else {
      segInt = k * std::log(E1/E0);
    }
    if (segInt < 0.0) segInt = 0.0;
    m_specCum[i+1] = m_specCum[i] + segInt;
  }

  // Normalize CDF
  const G4double total = m_specCum.back();
  if (total <= 0.0) {
    G4Exception("GPD3D_PrimaryGeneratorAction::LoadSpectrumFromCsv",
                "GPD3D_SPECTRUM_0003",
                JustWarning,
                "Spectrum integral is zero; check CSV weights.");
    m_specE.clear();
    m_specW.clear();
    m_specSegA.clear();
    m_specCum.clear();
    return;
  }
  for (auto& c : m_specCum) c /= total;
}

G4double GPD3D_PrimaryGeneratorAction::SampleSpectrumEnergy() const
{
  // Assumes m_specCum is normalized and non-empty.
  const G4double u = G4UniformRand();

  // Find segment i such that CDF[i] <= u < CDF[i+1]
  auto it = std::upper_bound(m_specCum.begin(), m_specCum.end(), u);
  std::size_t i = (it == m_specCum.begin()) ? 0 : (std::size_t)(it - m_specCum.begin() - 1);
  if (i >= m_specE.size() - 1) i = m_specE.size() - 2;

  const G4double c0 = m_specCum[i];
  const G4double c1 = m_specCum[i+1];
  // t is the local fraction of the segment: u=c0 means t=0, and u=c1 means t=1
  const G4double t = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0; // local [0,1]

  const G4double E0 = m_specE[i];
  const G4double E1 = m_specE[i+1];
  const G4double a  = m_specSegA[i];

  // Invert ∫E0..E W(E)dE for a power-law segment.
  if (std::fabs(a + 1.0) > 1e-14) {
    const G4double p0 = std::pow(E0, a + 1.0);
    const G4double p1 = std::pow(E1, a + 1.0);
    const G4double p  = p0 + t * (p1 - p0);
    return std::pow(p, 1.0 / (a + 1.0));
  }
  // a == -1 -> W(E) ~ 1/E, integral ~ ln(E)
  return E0 * std::exp(t * std::log(E1/E0));
}