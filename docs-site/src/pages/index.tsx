import Link from '@docusaurus/Link';
import Layout from '@theme/Layout';
import useBaseUrl from '@docusaurus/useBaseUrl';
import styles from './index.module.css';

const guideLinks = [
  {
    label: 'Get started',
    title: 'Install and run',
    text: 'Install from PyPI, Docker, or Apptainer and execute the shipped WES test.',
    to: '/docs/usage/quickstart',
  },
  {
    label: 'Analyses',
    title: 'Select a pipeline',
    text: 'Review the native WES, WGS, and mtDNA workflows for single-sample and cohort analysis.',
    to: '/docs/pipelines/overview',
  },
  {
    label: 'Reproducibility',
    title: 'Audit independent runs',
    text: 'Compare workflow identity, resources, runtime evidence, and final-output fingerprints.',
    to: '/docs/help/run-comparison',
  },
  {
    label: 'Developers',
    title: 'Add an implementation',
    text: 'Register a native pipeline or connect a supported external workflow provider.',
    to: '/docs/technical-details/extending-cbicall',
  },
];

export default function Home() {
  const logo = useBaseUrl('/img/cbicall-logo-transparent.png');
  const variantCallingVisual = useBaseUrl('/img/cbicall-variant-calling.svg');

  return (
    <Layout
      title="CBIcall"
      description="Configuration-driven execution and reproducibility auditing for variant-calling workflows">
      <main className={styles.page}>
        <section className={styles.hero}>
          <div className={styles.heroInner}>
            <div className={styles.copy}>
              <img className={styles.logo} src={logo} alt="CBIcall" />
              <h1 className={styles.srOnly}>CBIcall</h1>
              <p className={styles.kicker}>Configuration-driven variant calling</p>
              <p className={styles.lede}>
                Validate one analysis request against controlled workflow and
                resource contracts, execute it through supported backends, and
                retain structured evidence for audit and run comparison.
              </p>
              <div className={styles.actions}>
                <Link className="button button--primary button--lg" to="/docs/usage/quickstart">
                  Quick start
                </Link>
                <Link className="button button--secondary button--lg" to="/docs/overview">
                  Overview
                </Link>
              </div>
            </div>

            <div className={styles.contractPreview} aria-label="Example validated CBIcall analysis request">
              <div className={styles.previewTitle}>Analysis request</div>
              <pre><code><span>mode:</span> single{`\n`}<span>pipeline:</span> wes{`\n`}<span>workflow_backend:</span> nextflow{`\n`}<span>software_stack:</span> gatk-4.6{`\n`}<span>genome:</span> b37</code></pre>
              <div className={styles.resolution}>
                <div><span>Workflow</span><strong>resolved</strong></div>
                <div><span>Resources</span><strong>verified</strong></div>
                <div><span>Audit</span><strong>enabled</strong></div>
              </div>
            </div>
          </div>
        </section>

        <section className={styles.workflow} aria-label="Variant-calling workflow">
          <div className={styles.workflowInner}>
            <img
              className={styles.workflowImage}
              src={variantCallingVisual}
              alt="Paired-end FASTQ reads aligned into a BAM pileup and reported as a filtered VCF"
            />
            <div className={styles.mobileFlow}>
              <div><span>FASTQ</span><strong>Paired-end reads</strong></div>
              <div><span>BAM</span><strong>Aligned read evidence</strong></div>
              <div><span>VCF</span><strong>Filtered variant calls</strong></div>
            </div>
            <p>From paired-end sequencing reads to filtered VCFs with structured audit evidence.</p>
          </div>
        </section>

        <section className={styles.sections} aria-label="Documentation sections">
          <div className={styles.grid}>
            {guideLinks.map((guide) => (
              <Link className={styles.card} to={guide.to} key={guide.title}>
                <span>{guide.label}</span>
                <h2>{guide.title}</h2>
                <p>{guide.text}</p>
              </Link>
            ))}
          </div>
        </section>
      </main>
    </Layout>
  );
}
