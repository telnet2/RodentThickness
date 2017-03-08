#include <QObject>
#include <QProgress>

class PythonRunner : public QObject {
  Q_OBJECT

  public:
    PythonRunner(GuiCSV* parent) {
      _parent = parent;
    }

  public slots:
    void runScript() {
      if (_script.state() == QProcess::NotRunning) {
        QStringList args;
        args << __scriptPath;
        args << "--computeThickness";
        args << "--precorrespondence";
        args << "--runcorrespondence";
        args << "--compute-statistics";
        args
        t() << __scriptPath << "--computeThickness --precorrespondence --runcoorespondence --compute-statistics data.csv config.bms outputdir");
        _script.start(__pythonPath, args);
      }
    }
    void stopScript(QProcess::ProcessState state) {
      if (state == QProcess::NotRunning) {
        return;
      }
      _script.kill();
    }

  private:
    GuiCSV* _parent;
    QProcess _script;
}
