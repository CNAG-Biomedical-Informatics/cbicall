import Link from '@docusaurus/Link';
import Layout from '@theme/Layout';
import useBaseUrl from '@docusaurus/useBaseUrl';
import styles from './index.module.css';

export default function Home() {
  const logoUrl = useBaseUrl('/img/cbicall-logo-transparent.png');

  return (
    <Layout
      title="CBIcall"
      description="CNAG Biomedical Informatics framework for variant calling">
      <main>
        <section className={styles.hero}>
          <div className={styles.heroGrid}>
            <div className={styles.heroCopy}>
              <img className={styles.heroLogo} src={logoUrl} alt="CBIcall" />
              <h1 className={styles.srOnly}>CBIcall</h1>
              <p className={styles.heroSubtitle}>
                A configuration-driven workflow dispatcher for reproducible WES,
                WGS, and mtDNA variant-calling runs.
              </p>
              <div className={styles.actions}>
                <Link className="button button--primary button--lg" to="/docs/usage/quickstart">
                  Quick Start
                </Link>
                <Link className="button button--secondary button--lg" to="/docs/overview">
                  Overview
                </Link>
              </div>
            </div>

            <div className={styles.workflowPanel} aria-label="CBIcall workflow summary">
              <div className={styles.workflowStep}>
                <span>Input</span>
                <strong>FASTQ; existing BAMs for mtDNA</strong>
              </div>
              <div className={styles.workflowArrow}>→</div>
              <div className={styles.workflowStepPrimary}>
                <span>Dispatcher</span>
                <strong>WES · WGS · mtDNA</strong>
              </div>
              <div className={styles.workflowArrow}>→</div>
              <div className={styles.workflowStep}>
                <span>Output</span>
                <strong>BAM · VCF · QC · HTML</strong>
              </div>
              <div className={styles.tokens}>
                <span>bash</span>
                <span>snakemake</span>
                <span>nextflow</span>
                <span>GATK</span>
                <span>MToolBox</span>
              </div>
            </div>
          </div>
        </section>
        <section className={styles.sections}>
          <div className={styles.grid}>
            <Link className={styles.card} to="/docs/installation/docker">
              <span>Setup</span>
              <h2>Install</h2>
              <p>Docker, Apptainer, and non-containerized setup paths.</p>
            </Link>
            <Link className={styles.card} to="/docs/usage/quickstart">
              <span>Workflow</span>
              <h2>Run</h2>
              <p>Single-sample and cohort workflows for WES, WGS, and mtDNA.</p>
            </Link>
            <Link className={styles.card} to="/docs/help/outputs">
              <span>Results</span>
              <h2>Inspect</h2>
              <p>Standard project layout, logs, QC files, and generated outputs.</p>
            </Link>
            <Link className={styles.card} to="/docs/technical-details/adding-a-pipeline">
              <span>Developers</span>
              <h2>Extend</h2>
              <p>Registry-driven workflow wiring for adding new pipelines.</p>
            </Link>
          </div>
        </section>
      </main>
    </Layout>
  );
}
